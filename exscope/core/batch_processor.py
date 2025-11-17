"""Batch processing for multiple genes or regions.

This module provides functionality to process and visualize multiple genes
or genomic regions in parallel, utilizing multiprocessing for efficiency.
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional, Union
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from exscope.core.coverage import CoverageCalculator, SequencingType
from exscope.core.annotations import AnnotationProcessor
from exscope.visualizers.genomics_track import GenomeVisualizer, GenomicRegion


@dataclass
class GeneTarget:
    """Represents a gene target for visualization."""
    name: str
    region: Optional[GenomicRegion] = None
    
    def __str__(self) -> str:
        if self.region:
            return f"{self.name} ({self.region})"
        return self.name


class BatchProcessor:
    """Process and visualize multiple genes in batch mode."""
    
    def __init__(self, bam_path: Union[str, Path],
                 annotation_path: Union[str, Path],
                 output_dir: Union[str, Path],
                 reference: Optional[Union[str, Path]] = None,
                 threads: int = 4,
                 dpi: int = 100):
        """Initialize batch processor.
        
        Args:
            bam_path: Path to BAM/CRAM file
            annotation_path: Path to annotation file (BED/GTF/GFF)
            output_dir: Directory for output files
            reference: Reference genome (required for CRAM)
            threads: Number of parallel processes
            dpi: Image resolution
        """
        self.bam_path = Path(bam_path)
        self.annotation_path = Path(annotation_path)
        self.output_dir = Path(output_dir)
        self.reference = Path(reference) if reference else None
        self.threads = threads
        self.dpi = dpi
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize annotation processor once for all genes
        self.annotation_processor = AnnotationProcessor(self.annotation_path)
        
        logging.info(f"Batch processor initialized with {threads} threads")
    
    def extract_gene_regions(self, gene_names: List[str]) -> Dict[str, GenomicRegion]:
        """Extract genomic regions for multiple genes.
        
        Args:
            gene_names: List of gene names to extract
            
        Returns:
            Dictionary mapping gene names to their genomic regions
        """
        regions = {}
        
        for gene_name in gene_names:
            try:
                if hasattr(self.annotation_processor, 'db') and self.annotation_processor.db:
                    # GTF/GFF format
                    feature = self.annotation_processor.db[gene_name]
                    regions[gene_name] = GenomicRegion(
                        chrom=feature.seqid,
                        start=int(feature.start) if feature.start is not None else 0,
                        end=int(feature.end) if feature.end is not None else 0
                    )
                elif hasattr(self.annotation_processor, 'bedtool') and self.annotation_processor.bedtool:
                    # BED format - iterate through entries
                    filtered = []
                    for entry in self.annotation_processor.bedtool:
                        if entry.name == gene_name or entry.name.startswith(f"{gene_name}_"):
                            filtered.append(entry)
                    
                    if filtered:
                        regions[gene_name] = GenomicRegion(
                            chrom=filtered[0].chrom,
                            start=min(int(f.start) for f in filtered),
                            end=max(int(f.end) for f in filtered)
                        )
                    else:
                        logging.warning(f"Gene {gene_name} not found in annotation")
                else:
                    logging.warning(f"Unable to extract region for {gene_name}")
                    
            except Exception as e:
                logging.error(f"Failed to extract region for {gene_name}: {e}")
                
        return regions
    
    def _process_single_gene(self, gene_name: str, region: GenomicRegion,
                            coverage_prefix: str) -> Dict:
        """Process a single gene (internal method for parallel execution).
        
        Args:
            gene_name: Name of the gene
            region: Genomic region
            coverage_prefix: Prefix for coverage output files
            
        Returns:
            Dictionary with processing results
        """
        try:
            # Calculate coverage for this gene region
            calculator = CoverageCalculator(
                bam_path=self.bam_path,
                seq_type=SequencingType.PANEL,
                target_regions=self.annotation_path,
                reference=self.reference,
                threads=1  # Use 1 thread per gene since we're parallelizing genes
            )
            
            coverage_stats = calculator.calculate_coverage(coverage_prefix)
            
            # Get gene-specific coverage
            gene_coverage = calculator.get_gene_coverage(coverage_prefix, gene_name)
            
            # Create visualization
            visualizer = GenomeVisualizer(dpi=self.dpi)
            
            # Add coverage track
            coverage_file = Path(f"{coverage_prefix}.regions.bed.gz")
            if coverage_file.exists():
                visualizer.add_coverage_track(
                    coverage_file,
                    title=f"{gene_name} Coverage ({gene_coverage['coverage'].mean():.1f}x)"
                )
            
            # Add annotation track
            visualizer.add_coverage_track(
                self.annotation_path,
                title="Target Regions",
                color="#1F78B4"
            )
            
            # Generate visualization
            output_file = self.output_dir / f"{gene_name}_visualization.png"
            visualizer.plot(region=region, output=output_file, width=40, height_per_track=2)
            
            return {
                'gene': gene_name,
                'status': 'success',
                'output': str(output_file),
                'mean_coverage': float(gene_coverage['coverage'].mean()),
                'median_coverage': float(gene_coverage['coverage'].median()),
                'regions': len(gene_coverage)
            }
            
        except Exception as e:
            logging.error(f"Failed to process {gene_name}: {e}")
            return {
                'gene': gene_name,
                'status': 'failed',
                'error': str(e)
            }
    
    def process_genes(self, gene_names: List[str], parallel: bool = True) -> List[Dict]:
        """Process multiple genes and generate visualizations.
        
        Args:
            gene_names: List of gene names to process
            parallel: Whether to use parallel processing
            
        Returns:
            List of processing results for each gene
        """
        logging.info(f"Processing {len(gene_names)} genes in {'parallel' if parallel else 'sequential'} mode")
        
        # Extract regions for all genes
        gene_regions = self.extract_gene_regions(gene_names)
        
        if not gene_regions:
            raise ValueError("No valid gene regions found")
        
        results = []
        
        if parallel and len(gene_regions) > 1:
            # Parallel processing
            max_workers = min(self.threads, len(gene_regions))
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {}
                
                for gene_name, region in gene_regions.items():
                    coverage_prefix = str(self.output_dir / f"{gene_name}_coverage")
                    future = executor.submit(
                        self._process_single_gene,
                        gene_name, region, coverage_prefix
                    )
                    futures[future] = gene_name
                
                for future in as_completed(futures):
                    gene_name = futures[future]
                    try:
                        result = future.result()
                        results.append(result)
                        if result['status'] == 'success':
                            logging.info(f"✓ Completed {gene_name}")
                        else:
                            logging.warning(f"✗ Failed {gene_name}")
                    except Exception as e:
                        logging.error(f"✗ Exception processing {gene_name}: {e}")
                        results.append({
                            'gene': gene_name,
                            'status': 'failed',
                            'error': str(e)
                        })
        else:
            # Sequential processing
            for gene_name, region in gene_regions.items():
                coverage_prefix = str(self.output_dir / f"{gene_name}_coverage")
                result = self._process_single_gene(gene_name, region, coverage_prefix)
                results.append(result)
                
        # Log summary
        successful = sum(1 for r in results if r['status'] == 'success')
        failed = len(results) - successful
        
        logging.info(f"Batch processing complete: {successful} successful, {failed} failed")
        
        return results
    
    def generate_summary_report(self, results: List[Dict], output_file: Optional[Path] = None):
        """Generate a summary report of batch processing results.
        
        Args:
            results: List of processing results
            output_file: Optional path for summary report
        """
        if not output_file:
            output_file = self.output_dir / "batch_summary.txt"
        
        with open(output_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("Batch Processing Summary\n")
            f.write("=" * 60 + "\n\n")
            
            successful = [r for r in results if r['status'] == 'success']
            failed = [r for r in results if r['status'] == 'failed']
            
            f.write(f"Total genes: {len(results)}\n")
            f.write(f"Successful: {len(successful)}\n")
            f.write(f"Failed: {len(failed)}\n\n")
            
            if successful:
                f.write("Successful Genes:\n")
                f.write("-" * 60 + "\n")
                for r in successful:
                    f.write(f"  {r['gene']}: "
                           f"{r['mean_coverage']:.1f}x mean coverage, "
                           f"{r['regions']} regions\n")
                f.write("\n")
            
            if failed:
                f.write("Failed Genes:\n")
                f.write("-" * 60 + "\n")
                for r in failed:
                    f.write(f"  {r['gene']}: {r.get('error', 'Unknown error')}\n")
        
        logging.info(f"Summary report saved to: {output_file}")
