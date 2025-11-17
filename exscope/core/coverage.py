"""Coverage calculation using mosdepth.

This module provides a wrapper around mosdepth for efficient coverage calculation
with support for different NGS data types (WGS, WES, targeted panels) and output formats.
"""

import subprocess
from pathlib import Path
import logging
import pandas as pd
from typing import Optional, Union, Dict, List
import numpy as np
from enum import Enum

class SequencingType(Enum):
    """Types of NGS sequencing."""
    WGS = "whole_genome"
    WES = "whole_exome"
    PANEL = "targeted_panel"

class CoverageCalculator:
    """Wrapper for mosdepth coverage calculation with NGS-specific optimizations."""
    
    def __init__(self, bam_path: Union[str, Path], seq_type: SequencingType,
                 target_regions: Optional[Union[str, Path]] = None,
                 reference: Optional[Union[str, Path]] = None,
                 threads: int = 4):
        """Initialize CoverageCalculator.
        
        Args:
            bam_path: Path to the BAM/CRAM file
            seq_type: Type of sequencing data (WGS, WES, or targeted panel)
            target_regions: BED file with target regions (required for WES/panel)
            reference: Optional path to reference genome (required for CRAM)
            threads: Number of threads to use
        """
        self.bam_path = Path(bam_path)
        self.seq_type = seq_type
        self.target_regions = Path(target_regions) if target_regions else None
        self.reference = Path(reference) if reference else None
        self.threads = threads
        self._validate_files()
        
    def _validate_files(self):
        """Validate input files and sequencing type requirements."""
        if not self.bam_path.exists():
            raise FileNotFoundError(f"BAM/CRAM file not found: {self.bam_path}")
            
        if self.bam_path.suffix.lower() == '.cram' and not self.reference:
            raise ValueError("Reference genome is required for CRAM files")
            
        if self.reference and not self.reference.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference}")
            
        # Validate target regions requirement
        if self.seq_type in [SequencingType.WES, SequencingType.PANEL]:
            if not self.target_regions:
                raise ValueError(f"Target regions BED file required for {self.seq_type.value}")
            if not self.target_regions.exists():
                raise FileNotFoundError(f"Target regions file not found: {self.target_regions}")
                
    def calculate_coverage(self, output_prefix: str, window_size: Optional[int] = None,
                         min_coverage: int = 20) -> Dict:
        """Calculate coverage using mosdepth with sequencing type-specific settings.
        
        Args:
            output_prefix: Prefix for output files
            window_size: Window size for coverage calculation (default: auto-selected)
            min_coverage: Minimum coverage threshold for statistics
            
        Returns:
            Dictionary with coverage statistics
        """
        cmd = ['mosdepth', '-n', '-t', str(self.threads)]
        
        # Add sequencing type specific parameters
        if self.seq_type == SequencingType.WGS:
            # For WGS, use larger windows for efficiency
            window_size = window_size or 500
            cmd.extend(['-w', str(window_size)])
        else:
            # For WES/Panel, use target regions and per-base coverage
            cmd.extend(['-b', str(self.target_regions)])
            
        if self.reference:
            cmd.extend(['-f', str(self.reference)])
            
        # Set thresholds for coverage distribution
        thresholds = "1,10,20,30,50,100"
        cmd.extend(['-T', thresholds])
        
        cmd.extend([output_prefix, str(self.bam_path)])
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return self._analyze_coverage(output_prefix, min_coverage)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to calculate coverage: {e.stderr.decode()}")
            
    def _analyze_coverage(self, output_prefix: str, min_coverage: int) -> Dict:
        """Analyze coverage results and generate statistics.
        
        Args:
            output_prefix: Prefix used for mosdepth output files
            min_coverage: Minimum coverage threshold
        """
        # Load coverage data
        coverage_data = self.load_coverage(output_prefix)
        
        # Calculate statistics
        stats = {
            'mean_coverage': coverage_data['coverage'].mean(),
            'median_coverage': coverage_data['coverage'].median(),
            'coverage_uniformity': self._calculate_uniformity(coverage_data['coverage']),
            'pct_above_min': (coverage_data['coverage'] >= min_coverage).mean() * 100
        }
        
        # Add sequencing type specific metrics
        if self.seq_type != SequencingType.WGS:
            stats.update(self._calculate_target_metrics(coverage_data))
            
        return stats
        
    def _calculate_uniformity(self, coverage: pd.Series) -> float:
        """Calculate coverage uniformity (% within ±20% of mean)."""
        mean_cov = coverage.mean()
        lower = mean_cov * 0.8
        upper = mean_cov * 1.2
        return ((coverage >= lower) & (coverage <= upper)).mean() * 100
        
    def _calculate_target_metrics(self, coverage_data: pd.DataFrame) -> Dict:
        """Calculate target-specific metrics for WES/Panel data."""
        total_bases = (coverage_data['end'] - coverage_data['start']).sum()
        covered_bases = ((coverage_data['end'] - coverage_data['start'])[coverage_data['coverage'] > 0]).sum()
        
        return {
            'target_bases': total_bases,
            'covered_bases': covered_bases,
            'target_coverage': (covered_bases / total_bases * 100) if total_bases > 0 else 0
        }
        
    def load_coverage(self, output_prefix: str) -> pd.DataFrame:
        """Load coverage results into pandas DataFrame.
        
        Args:
            output_prefix: Prefix used for mosdepth output files
            
        Returns:
            DataFrame with columns: chrom, start, end, name, coverage
        """
        coverage_file = Path(f"{output_prefix}.regions.bed.gz")
            
        if not coverage_file.exists():
            raise FileNotFoundError(f"Coverage file not found: {coverage_file}")
            
        # Read mosdepth output format: chrom start end name coverage
        df = pd.read_csv(
            coverage_file,
            sep='\t',
            names=['chrom', 'start', 'end', 'name', 'coverage'],
            compression='gzip',
            comment='#'  # Skip comment/header lines
        )
        
        # Ensure proper data types
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df['coverage'] = df['coverage'].astype(float)
        
        return df
    
    def get_gene_coverage(self, output_prefix: str, gene_name: str) -> pd.DataFrame:
        """Get coverage data for a specific gene.
        
        Args:
            output_prefix: Prefix used for mosdepth output files
            gene_name: Name of the gene to extract coverage for
            
        Returns:
            DataFrame with coverage data for the specified gene
        """
        df = self.load_coverage(output_prefix)
        
        # Filter for the gene (handles both exact matches and partial matches like BRCA1_1)
        gene_df = df[df['name'].str.contains(gene_name, regex=False)].copy()
        
        if gene_df.empty:
            raise ValueError(f"No coverage data found for gene: {gene_name}")
            
        return gene_df