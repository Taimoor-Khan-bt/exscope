#!/usr/bin/env python3
"""Example script demonstrating exscope functionalities.

This script showcases the core features of exscope using example data:
1. BAM processing with samtools
2. Coverage calculation with mosdepth
3. Annotation processing
4. Genomic visualization
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt

from exscope.core.coverage import CoverageCalculator, SequencingType
from exscope.core.annotations import AnnotationProcessor
from exscope.visualizers.genomics_track import (
    GenomeVisualizer, 
    TrackConfig,
    GenomicRegion
)

# Set up paths
BASE_DIR = Path("/home/taimoor/genomics/personal tools")
DATA_DIR = BASE_DIR / "ex-data"
OUTPUT_DIR = BASE_DIR / "example_output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Input files
BAM_FILE = DATA_DIR / "Example_deduped.bam"
TARGET_BED = BASE_DIR / "xgen-exome-hyb-panel-v2-targets-hg38.bed"

def test_bam_processing():
    """Test BAM processing functionality."""
    print("\n=== Testing BAM Processing ===")
    print("\nBAM processing will be implemented in a future version.")

def test_coverage_calculation():
    """Test coverage calculation functionality."""
    print("\n=== Testing Coverage Calculation ===")
    
    from exscope.core.coverage import SequencingType
    
    # Initialize coverage calculator
    print("\nInitializing coverage calculator...")
    cov_calc = CoverageCalculator(
        bam_path=BAM_FILE,
        seq_type=SequencingType.PANEL,
        target_regions=TARGET_BED
    )
    
    # Calculate coverage for target regions
    print("\nCalculating coverage for target regions...")
    coverage_prefix = str(OUTPUT_DIR / "coverage")
    target_stats = cov_calc.calculate_coverage(
        output_prefix=coverage_prefix,
        min_coverage=20
    )
    
    print("\nCoverage Statistics:")
    for metric, value in target_stats.items():
        print(f"{metric}: {value}")
    
    # Calculate coverage for BRCA1 region
    print("\nCalculating coverage for BRCA1 region...")
    brca1_prefix = str(OUTPUT_DIR / "brca1_coverage")
    brca1_stats = cov_calc.calculate_coverage(
        output_prefix=brca1_prefix,
        min_coverage=20
    )
    
    print("\nBRCA1 Coverage Statistics:")
    for metric, value in brca1_stats.items():
        print(f"{metric}: {value}")

def test_annotation_processing():
    """Test annotation processing functionality."""
    print("\n=== Testing Annotation Processing ===")
    
    # Initialize annotation processor
    print("\nProcessing target regions...")
    with AnnotationProcessor(annotation_path=TARGET_BED) as anno_proc:
        # Calculate basic statistics
        print("\nCalculating target region statistics...")
        if anno_proc.bedtool is None:
            raise RuntimeError("Failed to initialize bedtool")
            
        stats = {
            'total_regions': len(list(anno_proc.bedtool)),
            'total_bases': sum(int(feature.stop) - int(feature.start) 
                             for feature in anno_proc.bedtool)
        }
        
        print("\nTarget Region Statistics:")
        for key, value in stats.items():
            print(f"{key}: {value}")
        
        # Extract BRCA1 region
        print("\nExtracting BRCA1 region...")
        # Create temporary BED file for BRCA1 region
        brca1_region = OUTPUT_DIR / "brca1_region.bed"
        with open(brca1_region, 'w') as f:
            f.write("chr17\t43044295\t43170245\tBRCA1\n")
        
        # Find overlaps
        brca1_targets = anno_proc.find_overlaps(brca1_region)
        
        # Save BRCA1 targets
        brca1_bed = OUTPUT_DIR / "brca1_targets.bed"
        if brca1_targets:
            brca1_targets.saveas(str(brca1_bed))
            print(f"\nBRCA1 target regions saved to: {brca1_bed}")
            
        # Clean up temporary file
        brca1_region.unlink()

def test_visualization():
    """Test genomic visualization functionality."""
    print("\n=== Testing Visualization ===")
    
    # Create visualization for BRCA1 region
    region = GenomicRegion("chr17", 43044295, 43170245)
    
    # Initialize visualizer
    vis = GenomeVisualizer(dpi=100)
    
    # Add tracks
    print("\nCreating visualization tracks...")
    
    # Add tracks if files exist
    
    # Coverage track
    coverage_file = OUTPUT_DIR / "coverage.regions.bed.gz"
    if coverage_file.exists():
        vis.add_coverage_track(
            coverage_file,
            title="Coverage",
            color="#33A02C"
        )
    
    # Target regions track
    target_file = OUTPUT_DIR / "brca1_targets.bed"
    if target_file.exists():
        vis.add_track(TrackConfig(
            title="Target Regions",
            file_path=target_file,
            track_type="bed",
            color="#1F78B4"
        ))
    
    # BAM alignment track
    bam_file = DATA_DIR / "Example_deduped.bam"
    if bam_file.exists() and bam_file.with_suffix(bam_file.suffix + '.bai').exists():
        vis.add_alignment_track(
            bam_file,
            title="Alignments",
            height=5
        )
    
    if not vis.tracks:
        print("\nNo valid tracks found. Make sure the required files exist.")
    
    # Create visualization
    print("\nGenerating visualization...")
    output_file = OUTPUT_DIR / "brca1_visualization.png"
    vis.plot(
        region,
        output=output_file,
        width=15,
        height_per_track=3
    )
    
    print(f"\nVisualization saved to: {output_file}")

def main():
    """Run all tests."""
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"\nOutput directory: {OUTPUT_DIR}")
    
    # Run tests
    test_bam_processing()
    test_coverage_calculation()
    test_annotation_processing()
    test_visualization()

if __name__ == "__main__":
    main()