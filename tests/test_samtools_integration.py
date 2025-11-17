"""Test samtools integration with BamProcessor."""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from exscope.core.bam_processor import BamProcessor

def test_samtools_integration():
    """Test samtools BAM file processing."""
    
    # Setup paths
    project_root = Path(__file__).parent.parent.parent
    bam_path = project_root / "ex-data/Example_deduped.bam"
    
    print("=" * 60)
    print("Testing Samtools Integration")
    print("=" * 60)
    
    # Initialize processor
    print("\n1. Initializing BamProcessor...")
    try:
        processor = BamProcessor(bam_path=bam_path, threads=4)
        print("✓ BamProcessor initialized successfully")
    except Exception as e:
        print(f"✗ BamProcessor initialization failed: {e}")
        return False
    
    # Test stats
    print("\n2. Getting BAM file statistics...")
    try:
        stats = processor.get_stats()
        print("✓ Statistics retrieved successfully")
        print(f"\nKey Statistics:")
        for key in ['sequences', 'reads mapped', 'reads unmapped', 'average length']:
            if key in stats:
                print(f"  {key}: {stats[key]}")
    except Exception as e:
        print(f"✗ Getting statistics failed: {e}")
        return False
    
    # Test region extraction (BRCA1 on chr17)
    print("\n3. Testing region extraction (chr17:43000000-43100000)...")
    try:
        reads_output = processor.view_region("chr17:43000000-43100000", min_quality=20)
        if reads_output:
            num_reads = len([l for l in reads_output.split('\n') if l and not l.startswith('@')])
            print(f"✓ Extracted {num_reads} reads from region")
        else:
            print("✓ Region extraction completed (no reads found)")
    except Exception as e:
        print(f"✗ Region extraction failed: {e}")
        return False
    
    # Test filtered extraction
    print("\n4. Testing filtered extraction (exclude unmapped reads)...")
    try:
        reads_output = processor.view_region(
            "chr17:43000000-43100000",
            min_quality=20,
            exclude_flags="0x4"  # Exclude unmapped reads
        )
        if reads_output:
            num_reads = len([l for l in reads_output.split('\n') if l and not l.startswith('@')])
            print(f"✓ Extracted {num_reads} mapped reads from region")
        else:
            print("✓ Filtered extraction completed (no reads found)")
    except Exception as e:
        print(f"✗ Filtered extraction failed: {e}")
        return False
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
    return True

if __name__ == "__main__":
    success = test_samtools_integration()
    sys.exit(0 if success else 1)
