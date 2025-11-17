#!/usr/bin/env python3
"""Example script demonstrating exscope core functionalities."""

import os
from pathlib import Path

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
BED_FILE = BASE_DIR / "xgen-exome-hyb-panel-v2-targets-hg38.bed"

def create_brca1_visualization():
    """Create visualization for BRCA1 region."""
    print("\n=== Creating BRCA1 Region Visualization ===")
    
    # Define BRCA1 region
    region = GenomicRegion("chr17", 43044295, 43170245)  # hg38 coordinates
    
    # Initialize visualizer
    vis = GenomeVisualizer()
    
    # Add tracks
    print("\nAdding visualization tracks...")
    
    # BAM alignment track
    vis.add_alignment_track(
        BAM_FILE,
        title="Alignments"
    )
    
    # Target regions track
    vis.add_track(TrackConfig(
        title="Target Regions",
        file_path=BED_FILE,
        track_type="bed",
        color="#1F78B4"
    ))
    
    # Create visualization
    print("\nGenerating visualization...")
    output_file = OUTPUT_DIR / "brca1_visualization.png"
    vis.plot(region, output_file)
    
    print(f"\nVisualization saved to: {output_file}")

def main():
    """Run example."""
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    
    # Create visualization
    create_brca1_visualization()

if __name__ == "__main__":
    main()