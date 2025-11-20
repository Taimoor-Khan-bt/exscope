"""Chromosome-level coverage analysis for detecting deletions and aneuploidies.

This module provides tools to calculate and analyze coverage across entire chromosomes
to detect deletions (e.g., Y chromosome deletion), monosomies, and other chromosomal
abnormalities. Particularly useful for WES data where per-chromosome coverage comparison
is critical for detecting large structural variations.
"""

import subprocess
from pathlib import Path
import logging
import pandas as pd
from typing import Optional, Union, Dict, List, Tuple
import numpy as np
from dataclasses import dataclass
from enum import Enum


class ChromosomeStatus(Enum):
    """Status of chromosome coverage."""
    NORMAL = "normal"
    REDUCED = "reduced"
    DELETED = "deleted"
    ELEVATED = "elevated"
    PARTIAL = "partial"  # Partially sequenced (e.g., targeted regions only)
    NOT_SEQUENCED = "not_sequenced"  # Completely absent from BAM


class Karyotype(Enum):
    """Inferred sex chromosome karyotype."""
    XX = "XX (Female)"
    XY = "XY (Male)"
    X0 = "X0 (Turner syndrome)"
    XXY = "XXY (Klinefelter syndrome)"
    XYY = "XYY"
    XXX = "XXX"
    UNKNOWN = "Unknown"


@dataclass
class ChromosomeCoverage:
    """Coverage statistics for a single chromosome."""
    chrom: str
    mean_coverage: float
    median_coverage: float
    total_bases: int
    covered_bases: int
    normalized_ratio: float
    status: ChromosomeStatus
    percent_covered: float = 0.0  # Percentage of chromosome actually sequenced
    covered_regions: Optional[List[Tuple[int, int]]] = None  # List of (start, end) tuples for visualization
    
    def __post_init__(self):
        """Initialize default values."""
        if self.covered_regions is None:
            self.covered_regions = []
    
    def __repr__(self):
        num_regions = len(self.covered_regions) if self.covered_regions else 0
        return (f"ChromosomeCoverage({self.chrom}: {self.mean_coverage:.2f}x, "
                f"ratio={self.normalized_ratio:.3f}, status={self.status.value}, "
                f"covered={self.percent_covered:.1f}%, regions={num_regions})")


@dataclass
class ChromosomeAnalysisResult:
    """Complete chromosome coverage analysis result."""
    sample_id: str
    chromosomes: Dict[str, ChromosomeCoverage]
    autosomal_mean: float
    inferred_karyotype: Karyotype
    deletions: List[str]
    warnings: List[str]
    
    def get_summary_table(self) -> pd.DataFrame:
        """Generate summary table of all chromosomes."""
        data = []
        for chrom in sorted(self.chromosomes.keys(), key=self._sort_chromosomes):
            cov = self.chromosomes[chrom]
            data.append({
                'Chromosome': cov.chrom,
                'Mean_Coverage': f"{cov.mean_coverage:.2f}",
                'Normalized_Ratio': f"{cov.normalized_ratio:.3f}",
                'Status': cov.status.value,
                'Covered_Bases': cov.covered_bases,
                'Total_Bases': cov.total_bases
            })
        return pd.DataFrame(data)
    
    @staticmethod
    def _sort_chromosomes(chrom: str) -> Tuple:
        """Sort chromosomes naturally (chr1, chr2, ..., chr22, chrX, chrY, chrM)."""
        chrom_clean = chrom.replace('chr', '')
        if chrom_clean.isdigit():
            return (0, int(chrom_clean))
        elif chrom_clean == 'X':
            return (1, 0)
        elif chrom_clean == 'Y':
            return (1, 1)
        elif chrom_clean in ['M', 'MT']:
            return (1, 2)
        else:
            return (2, chrom_clean)


class ChromosomeCoverageAnalyzer:
    """Analyze coverage at chromosome level to detect deletions and abnormalities."""
    
    # Thresholds for classification
    DELETION_THRESHOLD = 0.15  # <15% of autosomal mean = deleted
    REDUCED_THRESHOLD = 0.65   # <65% of autosomal mean = reduced
    ELEVATED_THRESHOLD = 1.45  # >145% of autosomal mean = elevated
    
    # Expected ratios for sex chromosomes (normalized to autosomes)
    MALE_X_RATIO = 0.5    # Males have 1 X vs 2 autosomes per cell
    MALE_Y_RATIO = 0.5    # Males have 1 Y vs 2 autosomes per cell
    FEMALE_X_RATIO = 1.0  # Females have 2 X vs 2 autosomes per cell
    FEMALE_Y_RATIO = 0.0  # Females have no Y chromosome
    
    def __init__(self, bam_path: Union[str, Path], 
                 reference: Optional[Union[str, Path]] = None,
                 threads: int = 2):
        """Initialize ChromosomeCoverageAnalyzer.
        
        Args:
            bam_path: Path to the BAM/CRAM file
            reference: Optional path to reference genome (required for CRAM)
            threads: Number of threads to use (default: 2 to prevent system overload)
        """
        self.bam_path = Path(bam_path)
        self.reference = Path(reference) if reference else None
        # Limit threads to prevent system overload
        self.threads = min(threads, 2)
        self._validate_files()
        
    def _validate_files(self):
        """Validate input files."""
        if not self.bam_path.exists():
            raise FileNotFoundError(f"BAM/CRAM file not found: {self.bam_path}")
            
        if self.bam_path.suffix.lower() == '.cram' and not self.reference:
            raise ValueError("Reference genome is required for CRAM files")
            
        if self.reference and not self.reference.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference}")
    
    def calculate_chromosome_coverage(self, output_prefix: str,
                                     chromosomes: Optional[List[str]] = None) -> str:
        """Calculate per-chromosome coverage using pysam to get actual covered regions.
        
        This method reads the BAM to determine which genomic positions actually have
        coverage, allowing accurate visualization of exome/targeted sequencing gaps.
        
        Args:
            output_prefix: Prefix for output files
            chromosomes: Optional list of chromosomes to analyze (default: chr1-22, X, Y, M)
            
        Returns:
            Path to the summary file
        """
        # Default human chromosomes
        if chromosomes is None:
            chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        
        logging.info(f"Calculating chromosome coverage with region tracking")
        
        try:
            import pysam
            
            # Get expected chromosome lengths for humans (hg38)
            expected_lengths = {
                'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
                'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
                'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
                'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
                'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
                'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
                'chrM': 16569
            }
            
            bam = pysam.AlignmentFile(str(self.bam_path))
            summary_data = []
            covered_regions = {}  # Store covered regions for visualization
            
            for chrom in chromosomes:
                try:
                    # Get all mapped reads for this chromosome
                    reads = [r for r in bam.fetch(chrom) if not r.is_unmapped]
                    
                    if len(reads) == 0:
                        # No coverage - chromosome not sequenced
                        expected_length = expected_lengths.get(chrom, 0)
                        summary_data.append({
                            'chrom': chrom,
                            'length': expected_length,
                            'bases': 0,
                            'mean': 0.0,
                            'min': 0,
                            'max': 0,
                            'percent_covered': 0.0,
                            'covered_regions': []
                        })
                        covered_regions[chrom] = []
                        continue
                    
                    # Calculate actual covered bases using pileup
                    # More accurate but slower - counts bases with >0 coverage
                    covered_bases = set()
                    total_depth = 0
                    
                    for pileupcolumn in bam.pileup(chrom, truncate=True):
                        covered_bases.add(pileupcolumn.pos)
                        total_depth += pileupcolumn.n
                    
                    num_covered_bases = len(covered_bases)
                    expected_length = expected_lengths.get(chrom, max(covered_bases) if covered_bases else 0)
                    
                    # Calculate coverage percentage
                    percent_covered = (num_covered_bases / expected_length * 100) if expected_length > 0 else 0
                    mean_coverage = (total_depth / num_covered_bases) if num_covered_bases > 0 else 0
                    
                    # Get covered regions for visualization (merge nearby positions)
                    regions = self._merge_positions_to_regions(covered_bases, max_gap=1000)
                    
                    summary_data.append({
                        'chrom': chrom,
                        'length': expected_length,
                        'bases': num_covered_bases,
                        'mean': mean_coverage,
                        'min': 0,  # Would need separate calculation
                        'max': 0,  # Would need separate calculation
                        'percent_covered': percent_covered,
                        'covered_regions': regions
                    })
                    covered_regions[chrom] = regions
                    
                    logging.info(f"{chrom}: {num_covered_bases:,} bases covered ({percent_covered:.2f}%), "
                               f"{len(regions)} regions, mean depth {mean_coverage:.2f}x")
                    
                except Exception as e:
                    # Chromosome not in BAM
                    logging.debug(f"{chrom} not found in BAM: {e}")
                    expected_length = expected_lengths.get(chrom, 0)
                    summary_data.append({
                        'chrom': chrom,
                        'length': expected_length,
                        'bases': 0,
                        'mean': 0.0,
                        'min': 0,
                        'max': 0,
                        'percent_covered': 0.0,
                        'covered_regions': []
                    })
                    covered_regions[chrom] = []
            
            bam.close()
            
            # Write summary file
            summary_file = f"{output_prefix}.mosdepth.summary.txt"
            with open(summary_file, 'w') as f:
                f.write("chrom\tlength\tbases\tmean\tmin\tmax\tpercent_covered\n")
                for data in summary_data:
                    f.write(f"{data['chrom']}\t{data['length']}\t{data['bases']}\t"
                           f"{data['mean']:.2f}\t{data['min']}\t{data['max']}\t"
                           f"{data['percent_covered']:.2f}\n")
            
            # Write covered regions file for visualization
            regions_file = f"{output_prefix}.covered_regions.tsv"
            with open(regions_file, 'w') as f:
                f.write("chrom\tstart\tend\n")
                for chrom, regions in covered_regions.items():
                    for start, end in regions:
                        f.write(f"{chrom}\t{start}\t{end}\n")
            
            logging.info(f"Wrote coverage summary to: {summary_file}")
            logging.info(f"Wrote covered regions to: {regions_file}")
            
            return summary_file
            
        except Exception as e:
            logging.error(f"Failed to calculate chromosome coverage: {e}")
            raise RuntimeError(f"Failed to calculate chromosome coverage: {e}")
    
    def parse_chromosome_summary(self, summary_file: str) -> pd.DataFrame:
        """Parse mosdepth summary output to extract per-chromosome statistics.
        
        Args:
            summary_file: Path to mosdepth summary file
            
        Returns:
            DataFrame with columns: chrom, length, bases, mean, min, max
        """
        df = pd.read_csv(summary_file, sep='\t')
        
        # Filter for chromosome-level data (not total_region or total)
        df = df[df['chrom'] != 'total'].copy()
        df = df[df['chrom'] != 'total_region'].copy()
        
        logging.info(f"Loaded coverage for {len(df)} chromosomes")
        logging.debug(f"Chromosomes: {df['chrom'].tolist()}")
        
        return df
    
    def analyze_chromosomes(self, output_prefix: str, 
                          sample_id: Optional[str] = None) -> ChromosomeAnalysisResult:
        """Perform complete chromosome-level coverage analysis.
        
        Args:
            output_prefix: Prefix for mosdepth output files
            sample_id: Optional sample identifier
            
        Returns:
            ChromosomeAnalysisResult with complete analysis
        """
        if sample_id is None:
            sample_id = self.bam_path.stem
        
        logging.info(f"Analyzing chromosome coverage for sample: {sample_id}")
        
        # Calculate coverage
        summary_file = self.calculate_chromosome_coverage(output_prefix)
        
        # Load covered regions for visualization
        regions_file = f"{output_prefix}.covered_regions.tsv"
        regions_by_chrom = {}
        try:
            regions_df = pd.read_csv(regions_file, sep='\t')
            for chrom in regions_df['chrom'].unique():
                chrom_regions = regions_df[regions_df['chrom'] == chrom]
                regions_by_chrom[chrom] = list(zip(chrom_regions['start'], chrom_regions['end']))
        except Exception as e:
            logging.warning(f"Could not load covered regions file: {e}")
        
        # Parse results
        coverage_df = self.parse_chromosome_summary(summary_file)
        
        # Calculate autosomal mean (chr1-22)
        autosomes = coverage_df[coverage_df['chrom'].str.match(r'chr\d+$')]
        if len(autosomes) == 0:
            # Try without 'chr' prefix
            autosomes = coverage_df[coverage_df['chrom'].str.match(r'^\d+$')]
        
        if len(autosomes) == 0:
            raise ValueError("No autosomal chromosomes found in coverage data")
        
        autosomal_mean = autosomes['mean'].mean()
        logging.info(f"Autosomal mean coverage: {autosomal_mean:.2f}x")
        
        # Define all expected human chromosomes
        expected_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        
        # Analyze each chromosome
        chromosomes = {}
        warnings = []
        
        # First, process chromosomes that are in the data
        for _, row in coverage_df.iterrows():
            chrom = row['chrom']
            mean_cov = row['mean']
            percent_covered = row.get('percent_covered', 100.0)
            
            # Calculate normalized ratio
            normalized_ratio = mean_cov / autosomal_mean if autosomal_mean > 0 else 0
            
            # Determine status based on coverage ratio and percent sequenced
            if percent_covered < 50:  # Less than 50% of chromosome sequenced
                status = ChromosomeStatus.PARTIAL
                warnings.append(f"{chrom} partially sequenced ({percent_covered:.1f}% of chromosome)")
            else:
                # Full or near-full chromosome - classify by coverage depth
                status = self._classify_coverage(normalized_ratio)
                
                # Generate warnings for abnormal coverage
                if status == ChromosomeStatus.DELETED:
                    warnings.append(f"{chrom} appears deleted (coverage ratio: {normalized_ratio:.3f})")
                elif status == ChromosomeStatus.REDUCED:
                    warnings.append(f"{chrom} has reduced coverage (ratio: {normalized_ratio:.3f})")
                elif status == ChromosomeStatus.ELEVATED:
                    warnings.append(f"{chrom} has elevated coverage (ratio: {normalized_ratio:.3f})")
            
            chromosomes[chrom] = ChromosomeCoverage(
                chrom=chrom,
                mean_coverage=mean_cov,
                median_coverage=row.get('median', mean_cov),
                total_bases=int(row['length']),
                covered_bases=int(row['bases']),
                normalized_ratio=normalized_ratio,
                status=status,
                percent_covered=percent_covered,
                covered_regions=regions_by_chrom.get(chrom, [])
            )
        
        # Add missing chromosomes as "not sequenced" (show in visualization)
        for expected_chrom in expected_chromosomes:
            if expected_chrom not in chromosomes:
                chromosomes[expected_chrom] = ChromosomeCoverage(
                    chrom=expected_chrom,
                    mean_coverage=0.0,
                    median_coverage=0.0,
                    total_bases=0,
                    covered_bases=0,
                    normalized_ratio=0.0,
                    status=ChromosomeStatus.NOT_SEQUENCED,
                    percent_covered=0.0,
                    covered_regions=[]
                )
                warnings.append(f"{expected_chrom} not sequenced (no data in BAM)")
        
        # Infer karyotype and detect deletions
        karyotype = self._infer_karyotype(chromosomes)
        deletions = [chrom for chrom, cov in chromosomes.items() 
                    if cov.status == ChromosomeStatus.DELETED and cov.mean_coverage > 0]
        
        logging.info(f"Inferred karyotype: {karyotype.value}")
        if deletions:
            logging.warning(f"Detected deletions: {', '.join(deletions)}")
        
        return ChromosomeAnalysisResult(
            sample_id=sample_id,
            chromosomes=chromosomes,
            autosomal_mean=autosomal_mean,
            inferred_karyotype=karyotype,
            deletions=deletions,
            warnings=warnings
        )
    
    def _merge_positions_to_regions(self, positions: set, max_gap: int = 1000) -> List[Tuple[int, int]]:
        """Merge covered positions into contiguous regions.
        
        Args:
            positions: Set of covered genomic positions
            max_gap: Maximum gap size to merge (default 1000bp)
            
        Returns:
            List of (start, end) tuples representing covered regions
        """
        if not positions:
            return []
        
        # Sort positions
        sorted_pos = sorted(positions)
        regions = []
        
        start = sorted_pos[0]
        prev = sorted_pos[0]
        
        for pos in sorted_pos[1:]:
            if pos - prev > max_gap:
                # Gap too large - start new region
                regions.append((start, prev + 1))
                start = pos
            prev = pos
        
        # Add final region
        regions.append((start, prev + 1))
        
        return regions
    
    def _classify_coverage(self, normalized_ratio: float) -> ChromosomeStatus:
        """Classify chromosome coverage status based on normalized ratio.
        
        Args:
            normalized_ratio: Coverage ratio relative to autosomal mean
            
        Returns:
            ChromosomeStatus classification
        """
        if normalized_ratio < self.DELETION_THRESHOLD:
            return ChromosomeStatus.DELETED
        elif normalized_ratio < self.REDUCED_THRESHOLD:
            return ChromosomeStatus.REDUCED
        elif normalized_ratio > self.ELEVATED_THRESHOLD:
            return ChromosomeStatus.ELEVATED
        else:
            return ChromosomeStatus.NORMAL
    
    def _infer_karyotype(self, chromosomes: Dict[str, ChromosomeCoverage]) -> Karyotype:
        """Infer sex chromosome karyotype from coverage ratios.
        
        Args:
            chromosomes: Dictionary of chromosome coverage data
            
        Returns:
            Inferred Karyotype
        """
        # Find X and Y chromosomes
        x_chrom = None
        y_chrom = None
        
        for chrom_name, cov in chromosomes.items():
            chrom_clean = chrom_name.replace('chr', '').upper()
            if chrom_clean == 'X':
                x_chrom = cov
            elif chrom_clean == 'Y':
                y_chrom = cov
        
        if not x_chrom:
            return Karyotype.UNKNOWN
        
        x_ratio = x_chrom.normalized_ratio
        y_ratio = y_chrom.normalized_ratio if y_chrom else 0.0
        
        logging.debug(f"X ratio: {x_ratio:.3f}, Y ratio: {y_ratio:.3f}")
        
        # Classify based on ratios
        # Y deleted or absent
        if y_ratio < self.DELETION_THRESHOLD:
            if 0.85 <= x_ratio <= 1.15:  # Normal X coverage
                return Karyotype.XX
            elif 0.4 <= x_ratio <= 0.6:  # Half X coverage
                return Karyotype.X0
            elif x_ratio > 1.4:  # Elevated X
                return Karyotype.XXX
        
        # Y present
        else:
            if 0.4 <= x_ratio <= 0.6 and 0.4 <= y_ratio <= 0.6:  # Both half
                return Karyotype.XY
            elif 0.85 <= x_ratio <= 1.15 and 0.4 <= y_ratio <= 0.6:  # X normal, Y half
                return Karyotype.XXY
            elif 0.4 <= x_ratio <= 0.6 and y_ratio > 0.8:  # X half, Y normal
                return Karyotype.XYY
        
        return Karyotype.UNKNOWN
    
    def compare_samples(self, results: List[ChromosomeAnalysisResult]) -> pd.DataFrame:
        """Compare chromosome coverage across multiple samples.
        
        Args:
            results: List of ChromosomeAnalysisResult objects
            
        Returns:
            DataFrame with comparison of normalized ratios across samples
        """
        if not results:
            raise ValueError("No results provided for comparison")
        
        # Get all chromosomes
        all_chroms = set()
        for result in results:
            all_chroms.update(result.chromosomes.keys())
        
        # Build comparison table
        data = []
        for chrom in sorted(all_chroms, key=ChromosomeAnalysisResult._sort_chromosomes):
            row = {'Chromosome': chrom}
            for result in results:
                if chrom in result.chromosomes:
                    cov = result.chromosomes[chrom]
                    row[f"{result.sample_id}_ratio"] = cov.normalized_ratio
                    row[f"{result.sample_id}_status"] = cov.status.value
                else:
                    row[f"{result.sample_id}_ratio"] = 0.0
                    row[f"{result.sample_id}_status"] = "missing"
            data.append(row)
        
        return pd.DataFrame(data)
