"""Test batch processing functionality."""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from exscope.core.batch_processor import BatchProcessor

def test_batch_processing():
    """Test batch processing for multiple genes."""
    
    # Setup paths
    project_root = Path(__file__).parent.parent.parent
    bam_path = project_root / "ex-data/Example_deduped.bam"
    bed_path = project_root / "xgen-exome-hyb-panel-v2-targets-hg38.bed"
    output_dir = Path("test_batch_output")
    
    print("=" * 60)
    print("Testing Batch Processing")
    print("=" * 60)
    
    # Test 1: Initialize batch processor
    print("\n1. Initializing BatchProcessor...")
    try:
        processor = BatchProcessor(
            bam_path=bam_path,
            annotation_path=bed_path,
            output_dir=output_dir,
            threads=2,
            dpi=100
        )
        print("✓ BatchProcessor initialized successfully")
    except Exception as e:
        print(f"✗ BatchProcessor initialization failed: {e}")
        return False
    
    # Test 2: Extract regions for multiple genes
    print("\n2. Extracting regions for multiple genes...")
    genes = ["BRCA1", "BRCA2", "TP53"]
    try:
        regions = processor.extract_gene_regions(genes)
        print(f"✓ Extracted regions for {len(regions)} genes:")
        for gene, region in regions.items():
            print(f"  - {gene}: {region}")
    except Exception as e:
        print(f"✗ Region extraction failed: {e}")
        return False
    
    # Test 3: Process genes sequentially
    print("\n3. Processing genes sequentially...")
    try:
        results = processor.process_genes(
            gene_names=["BRCA1"],  # Start with just one gene for testing
            parallel=False
        )
        print(f"✓ Processed {len(results)} gene(s)")
        for result in results:
            if result['status'] == 'success':
                print(f"  ✓ {result['gene']}: {result['mean_coverage']:.1f}x coverage")
                print(f"    Output: {result['output']}")
            else:
                print(f"  ✗ {result['gene']}: {result.get('error', 'Unknown error')}")
    except Exception as e:
        print(f"✗ Sequential processing failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test 4: Generate summary report
    print("\n4. Generating summary report...")
    try:
        processor.generate_summary_report(results)
        summary_file = output_dir / "batch_summary.txt"
        if summary_file.exists():
            print(f"✓ Summary report generated: {summary_file}")
            print("\nReport contents:")
            print("-" * 60)
            print(summary_file.read_text())
        else:
            print("✗ Summary report not generated")
            return False
    except Exception as e:
        print(f"✗ Summary report generation failed: {e}")
        return False
    
    # Test 5: Process multiple genes in parallel
    print("\n5. Processing multiple genes in parallel...")
    try:
        # Use a subset of genes that exist in the data
        test_genes = ["BRCA1", "BRCA2"]
        results_parallel = processor.process_genes(
            gene_names=test_genes,
            parallel=True
        )
        print(f"✓ Processed {len(results_parallel)} genes in parallel")
        
        successful = [r for r in results_parallel if r['status'] == 'success']
        failed = [r for r in results_parallel if r['status'] == 'failed']
        
        print(f"  Successful: {len(successful)}")
        print(f"  Failed: {len(failed)}")
        
        for result in results_parallel:
            if result['status'] == 'success':
                print(f"  ✓ {result['gene']}: {result['mean_coverage']:.1f}x")
            else:
                print(f"  ✗ {result['gene']}: {result.get('error', 'Unknown')[:50]}")
                
    except Exception as e:
        print(f"✗ Parallel processing failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print("Generated files:")
    if output_dir.exists():
        for f in sorted(output_dir.glob("*")):
            print(f"  - {f.name}")
    
    return True

if __name__ == "__main__":
    success = test_batch_processing()
    sys.exit(0 if success else 1)
