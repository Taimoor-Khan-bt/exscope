"""Tests for annotation processing functionality."""

import pytest
from pathlib import Path
import pybedtools
from exscope.core.annotations import AnnotationProcessor, AnnotationFormat, GenomicFeature

@pytest.fixture
def sample_gtf(tmp_path):
    """Create a sample GTF file."""
    gtf_content = [
        '##description: test GTF',
        'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "gene1"; transcript_id "trans1"',
        'chr1\ttest\tCDS\t150\t200\t.\t+\t.\tgene_id "gene1"; transcript_id "trans1"'
    ]
    gtf_path = tmp_path / "test.gtf"
    gtf_path.write_text('\n'.join(gtf_content))
    return gtf_path

@pytest.fixture
def sample_bed(tmp_path):
    """Create a sample BED file."""
    bed_content = [
        'chr1\t50\t150\tfeature1\t0\t+',
        'chr1\t175\t225\tfeature2\t0\t+'
    ]
    bed_path = tmp_path / "test.bed"
    bed_path.write_text('\n'.join(bed_content))
    return bed_path

def test_annotation_format_detection(sample_gtf, sample_bed):
    """Test file format detection."""
    gtf_proc = AnnotationProcessor(sample_gtf)
    assert gtf_proc.format == AnnotationFormat.GTF
    
    bed_proc = AnnotationProcessor(sample_bed)
    assert bed_proc.format == AnnotationFormat.BED

def test_gff_database_creation(sample_gtf):
    """Test GFF database creation and query."""
    proc = AnnotationProcessor(sample_gtf)
    proc.create_gff_db(force=True)
    
    # Check if database was created
    assert proc.db is not None
    
    # Query features
    exons = list(proc.db.features_of_type('exon'))
    assert len(exons) == 1
    assert exons[0].start == 100
    assert exons[0].end == 200

def test_feature_extraction(sample_gtf):
    """Test feature extraction with flanking regions."""
    proc = AnnotationProcessor(sample_gtf)
    output = Path("test_output.bed")
    
    try:
        # Extract exons with 50bp flanking
        proc.extract_features(
            feature_type=GenomicFeature.EXON,
            output=output,
            upstream=50,
            downstream=50
        )
        
        # Verify output
        assert output.exists()
        content = output.read_text()
        assert 'chr1' in content
        assert '50\t250' in content  # Original 100-200 + 50bp flanking
    finally:
        if output.exists():
            output.unlink()

def test_overlap_detection(sample_gtf, sample_bed):
    """Test overlap detection between annotation files."""
    proc = AnnotationProcessor(sample_gtf)
    
    # Find overlaps
    overlaps = proc.find_overlaps(
        sample_bed,
        min_overlap=0.1,
        strand_specific=True
    )
    
    # Convert to list for checking
    overlap_regions = list(overlaps)
    assert len(overlap_regions) > 0  # Should find at least one overlap
    
    # Check if overlaps are correct
    assert any(100 <= int(region.start) <= 200 for region in overlap_regions)

def test_feature_merging(sample_bed):
    """Test merging of nearby features."""
    proc = AnnotationProcessor(sample_bed)
    
    # Merge features within 25bp
    merged = proc.merge_features(distance=25)
    merged_regions = list(merged)
    
    # Should merge the two features (50-150 and 175-225)
    assert len(merged_regions) == 1  # Features should be merged
    assert int(merged_regions[0].start) == 50
    assert int(merged_regions[0].end) == 225

def test_invalid_file_handling():
    """Test handling of invalid files."""
    with pytest.raises(FileNotFoundError):
        AnnotationProcessor("nonexistent.gtf")
        
    # Test empty file
    empty_file = Path("empty.gtf")
    empty_file.touch()
    try:
        with pytest.raises(ValueError):
            AnnotationProcessor(empty_file)
    finally:
        empty_file.unlink()

def test_format_conversion(sample_gtf, tmp_path):
    """Test format conversion capabilities."""
    proc = AnnotationProcessor(sample_gtf)
    output = tmp_path / "converted.bed"
    
    # Convert GTF to BED
    proc.to_format(output, AnnotationFormat.BED)
    
    # Verify conversion
    assert output.exists()
    content = output.read_text().splitlines()
    assert len(content) > 0
    # Verify BED format (tab-separated, at least 6 columns)
    assert all(len(line.split('\t')) >= 6 for line in content)