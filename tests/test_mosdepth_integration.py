"""Test mosdepth integration with CoverageCalculator."""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from exscope.core.coverage import CoverageCalculator, SequencingType

def test_mosdepth_integration():
    """Test mosdepth coverage calculation."""
    
    # Setup paths (absolute paths from project root)
    project_root = Path(__file__).parent.parent.parent
    bam_path = project_root / "ex-data/Example_deduped.bam"
    bed_path = project_root / "xgen-exome-hyb-panel-v2-targets-hg38.bed"
    output_prefix = "test_coverage"
    
    print("=" * 60)
    print("Testing Mosdepth Integration")
    print("=" * 60)
    
    # Initialize calculator
    print("\n1. Initializing CoverageCalculator...")
    calculator = CoverageCalculator(
        bam_path=bam_path,
        seq_type=SequencingType.PANEL,
        target_regions=bed_path,
        threads=4
    )
    print("✓ CoverageCalculator initialized successfully")
    
    # Calculate coverage
    print("\n2. Calculating coverage with mosdepth...")
    try:
        stats = calculator.calculate_coverage(output_prefix, min_coverage=20)
        print("✓ Coverage calculation completed")
        print(f"\nCoverage Statistics:")
        for key, value in stats.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.2f}")
            else:
                print(f"  {key}: {value}")
    except Exception as e:
        print(f"✗ Coverage calculation failed: {e}")
        return False
    
    # Test loading coverage data
    print("\n3. Loading coverage data...")
    try:
        coverage_df = calculator.load_coverage(output_prefix)
        print(f"✓ Loaded {len(coverage_df)} regions")
        print(f"\nFirst 5 regions:")
        print(coverage_df.head())
    except Exception as e:
        print(f"✗ Loading coverage failed: {e}")
        return False
    
    # Test gene-specific coverage
    print("\n4. Testing gene-specific coverage extraction...")
    try:
        brca1_coverage = calculator.get_gene_coverage(output_prefix, "BRCA1")
        print(f"✓ Found {len(brca1_coverage)} regions for BRCA1")
        print(f"\nBRCA1 Coverage Summary:")
        print(f"  Mean coverage: {brca1_coverage['coverage'].mean():.2f}x")
        print(f"  Median coverage: {brca1_coverage['coverage'].median():.2f}x")
        print(f"  Max coverage: {brca1_coverage['coverage'].max():.2f}x")
        print(f"  Min coverage: {brca1_coverage['coverage'].min():.2f}x")
    except Exception as e:
        print(f"✗ Gene-specific coverage extraction failed: {e}")
        return False
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
    return True

if __name__ == "__main__":
    success = test_mosdepth_integration()
    sys.exit(0 if success else 1)
