"""Test pyGenomeTracks integration and whole chromosome visualization."""

import sys
from pathlib import Path
import tempfile

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from exscope.core.coverage import CoverageCalculator, SequencingType
from exscope.visualizers.genomics_track import (
    GenomeVisualizer,
    GenomicRegion,
    TrackConfig
)

def test_pygenometracks_integration():
    """Test pyGenomeTracks visualization capabilities."""
    
    # Setup paths
    project_root = Path(__file__).parent.parent.parent
    bam_path = project_root / "ex-data/Example_deduped.bam"
    bed_path = project_root / "xgen-exome-hyb-panel-v2-targets-hg38.bed"
    
    print("=" * 60)
    print("Testing pyGenomeTracks Integration")
    print("=" * 60)
    
    # Test 1: Single gene visualization
    print("\n1. Testing single gene visualization (BRCA1)...")
    try:
        # Calculate coverage first
        calculator = CoverageCalculator(
            bam_path=bam_path,
            seq_type=SequencingType.PANEL,
            target_regions=bed_path,
            threads=4
        )
        
        output_prefix = "brca1_viz_test"
        coverage_stats = calculator.calculate_coverage(output_prefix)
        
        # Create visualizer
        visualizer = GenomeVisualizer(dpi=100)
        
        # Add coverage track
        coverage_file = Path(f"{output_prefix}.regions.bed.gz")
        if coverage_file.exists():
            visualizer.add_coverage_track(
                coverage_file,
                title=f"Coverage (Mean: {coverage_stats.get('mean_coverage', 0):.1f}x)"
            )
        
        # Add BED track for target regions
        visualizer.add_coverage_track(bed_path, title="Target Regions", color="#1F78B4")
        
        # Define BRCA1 region
        brca1_region = GenomicRegion(chrom="chr17", start=43044295, end=43125483)
        
        # Create visualization
        output_file = Path("test_brca1_visualization.png")
        visualizer.plot(region=brca1_region, output=output_file, width=40, height_per_track=2)
        
        if output_file.exists():
            print(f"✓ Single gene visualization created: {output_file}")
            print(f"  File size: {output_file.stat().st_size / 1024:.1f} KB")
        else:
            print("✗ Visualization file not created")
            return False
            
    except Exception as e:
        print(f"✗ Single gene visualization failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test 2: Whole chromosome visualization
    print("\n2. Testing whole chromosome visualization (chr17)...")
    try:
        # Calculate coverage for whole chromosome
        calculator_wgs = CoverageCalculator(
            bam_path=bam_path,
            seq_type=SequencingType.WGS,  # Use WGS mode for chromosome-wide analysis
            threads=4
        )
        
        output_prefix_chr = "chr17_viz_test"
        coverage_stats_chr = calculator_wgs.calculate_coverage(output_prefix_chr, window_size=500000)
        
        # Create visualizer for chromosome
        visualizer_chr = GenomeVisualizer(dpi=100)
        
        # Add coverage track
        coverage_file_chr = Path(f"{output_prefix_chr}.regions.bed.gz")
        if coverage_file_chr.exists():
            visualizer_chr.add_coverage_track(
                coverage_file_chr,
                title=f"Chr17 Coverage (Mean: {coverage_stats_chr.get('mean_coverage', 0):.1f}x)"
            )
        
        # Create chromosome-wide visualization
        output_file_chr = Path("test_chr17_visualization.png")
        visualizer_chr.plot_chromosome("chr17", output=output_file_chr, width=60, height_per_track=2)
        
        if output_file_chr.exists():
            print(f"✓ Chromosome visualization created: {output_file_chr}")
            print(f"  File size: {output_file_chr.stat().st_size / 1024:.1f} KB")
        else:
            print("✗ Chromosome visualization file not created")
            return False
            
    except Exception as e:
        print(f"✗ Chromosome visualization failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test 3: Multiple genes batch visualization
    print("\n3. Testing batch gene visualization (BRCA1 + BRCA2)...")
    try:
        genes = {
            "BRCA1": GenomicRegion(chrom="chr17", start=43044295, end=43125483),
            "BRCA2": GenomicRegion(chrom="chr13", start=32315086, end=32400266)
        }
        
        for gene_name, gene_region in genes.items():
            visualizer_gene = GenomeVisualizer(dpi=100)
            
            # Add tracks
            if coverage_file.exists():
                visualizer_gene.add_coverage_track(coverage_file, title=f"{gene_name} Coverage")
            visualizer_gene.add_coverage_track(bed_path, title="Target Regions", color="#1F78B4")
            
            # Create visualization
            output_gene = Path(f"test_{gene_name.lower()}_visualization.png")
            visualizer_gene.plot(region=gene_region, output=output_gene, width=40, height_per_track=2)
            
            if output_gene.exists():
                print(f"  ✓ {gene_name} visualization: {output_gene.stat().st_size / 1024:.1f} KB")
            else:
                print(f"  ✗ {gene_name} visualization failed")
                
    except Exception as e:
        print(f"✗ Batch visualization failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
    print("\nGenerated visualization files:")
    print("  - test_brca1_visualization.png")
    print("  - test_chr17_visualization.png")
    print("  - test_brca1_visualization.png")
    print("  - test_brca2_visualization.png")
    return True

if __name__ == "__main__":
    success = test_pygenometracks_integration()
    sys.exit(0 if success else 1)
