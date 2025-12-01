"""Tests for genomic visualization."""

import pytest
from pathlib import Path
import matplotlib.pyplot as plt
import yaml

from exscope.visualizers.genomics_track import (
    GenomicRegion,
    TrackConfig,
    GenomeVisualizer
)

def test_genomic_region():
    """Test GenomicRegion class."""
    region = GenomicRegion("chr1", 1000, 2000)
    assert str(region) == "chr1:1000-2000"
    assert region.chrom == "chr1"
    assert region.start == 1000
    assert region.end == 2000

def test_track_config_validation():
    """Test TrackConfig validation."""
    # Test valid configurations
    track = TrackConfig(
        title="Coverage",
        file_path="test.bw",
        track_type="bigwig"
    )
    assert track.title == "Coverage"
    assert track.track_type == "bigwig"
    
    # Test invalid track type
    with pytest.raises(ValueError):
        TrackConfig(
            title="Invalid",
            file_path="test.txt",
            track_type="invalid"
        )
    
    # Test invalid file format
    with pytest.raises(ValueError):
        TrackConfig(
            title="Wrong format",
            file_path="test.txt",
            track_type="bigwig"
        )

def test_track_config_to_dict():
    """Test TrackConfig conversion to dict."""
    track = TrackConfig(
        title="Coverage",
        file_path="test.bw",
        track_type="bigwig",
        height=5,
        color="#FF0000"
    )
    config = track.to_dict()
    
    assert config["title"] == "Coverage"
    assert config["file"] == "test.bw"
    assert config["height"] == 5
    assert config["color"] == "#FF0000"
    assert config["file_type"] == "bigwig"

def test_genome_visualizer():
    """Test GenomeVisualizer class."""
    vis = GenomeVisualizer(dpi=100)
    
    # Add tracks
    vis.add_coverage_track("coverage.bw", title="Coverage")
    vis.add_annotation_track("genes.gtf", title="Genes")
    vis.add_alignment_track("reads.bam", title="Alignments")
    
    assert len(vis.tracks) == 3
    assert vis.tracks[0].track_type == "bigwig"
    assert vis.tracks[1].track_type == "gtf"
    assert vis.tracks[2].track_type == "bam"

def test_config_file_generation(tmp_path):
    """Test track configuration file generation."""
    vis = GenomeVisualizer()
    vis.add_coverage_track("coverage.bw", title="Coverage")
    
    config_file = tmp_path / "tracks.ini"
    vis.save_config(config_file)
    
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    assert "[track1]" in config
    assert config["[track1]"]["title"] == "Coverage"
    assert config["[track1]"]["file_type"] == "bigwig"

def test_plot_validation():
    """Test plot validation."""
    vis = GenomeVisualizer()
    region = GenomicRegion("chr1", 1000, 2000)
    
    # Test plotting without tracks
    with pytest.raises(ValueError):
        vis.plot(region)
        
    # Add track and test plotting
    vis.add_coverage_track("coverage.bw")
    try:
        vis.plot(region)
    except Exception as e:
        # Only file not found errors are acceptable
        assert "No such file" in str(e)


# ============================================================================
# Tests for Chromosome Coverage Visualization (v1.3.0+)
# ============================================================================

def test_get_ratio_color():
    """Test biological interpretation-based color assignment."""
    from exscope.visualizers.chromosome_plot import ChromosomePlotter
    
    plotter = ChromosomePlotter()
    
    # Test normal diploid (green)
    assert plotter.get_ratio_color(1.0) == '#2ecc71'
    assert plotter.get_ratio_color(0.95) == '#2ecc71'
    assert plotter.get_ratio_color(1.05) == '#2ecc71'
    
    # Test haploid (yellow/orange)
    assert plotter.get_ratio_color(0.5) == '#f39c12'
    assert plotter.get_ratio_color(0.45) == '#f39c12'
    assert plotter.get_ratio_color(0.55) == '#f39c12'
    
    # Test partial loss/gain (orange)
    assert plotter.get_ratio_color(0.75) == '#e67e22'
    assert plotter.get_ratio_color(1.2) == '#e67e22'
    
    # Test deletion/duplication (red)
    assert plotter.get_ratio_color(0.2) == '#e74c3c'
    assert plotter.get_ratio_color(1.6) == '#e74c3c'
    
    # Test not sequenced (gray)
    assert plotter.get_ratio_color(0.0) == '#95a5a6'


def test_get_ratio_interpretation():
    """Test clinical interpretation text generation."""
    from exscope.visualizers.chromosome_plot import ChromosomePlotter
    from exscope.core.chromosome_coverage import Karyotype
    
    plotter = ChromosomePlotter()
    
    # Test autosomal chromosomes
    assert "Normal diploid" in plotter.get_ratio_interpretation(1.0, "chr1", Karyotype.XX)
    assert "Hemizygous deletion" in plotter.get_ratio_interpretation(0.5, "chr2", Karyotype.XY)
    assert "Complete deletion" in plotter.get_ratio_interpretation(0.1, "chr3", Karyotype.XX)
    assert "Trisomy" in plotter.get_ratio_interpretation(1.5, "chr21", Karyotype.XY)
    assert "Tetrasomy" in plotter.get_ratio_interpretation(2.0, "chr4", Karyotype.XX)
    
    # Test sex chromosomes in males (XY)
    assert "Normal" in plotter.get_ratio_interpretation(0.5, "chrX", Karyotype.XY)
    assert "Normal" in plotter.get_ratio_interpretation(0.5, "chrY", Karyotype.XY)
    # For XXY, ratio would be 1.0 for chrX but karyotype context matters - skip this complex case
    
    # Test sex chromosomes in females (XX)
    assert "Normal" in plotter.get_ratio_interpretation(1.0, "chrX", Karyotype.XX)
    # chrY present in female is abnormal - but get_ratio_interpretation doesn't flag this, just shows ratio
    
    # Test Turner syndrome (X0) - chrX ratio of 0.5 shows as hemizygous deletion
    assert "Hemizygous" in plotter.get_ratio_interpretation(0.5, "chrX", Karyotype.X0)


def test_chromosome_analysis_result_creation():
    """Test creation of synthetic ChromosomeAnalysisResult for testing."""
    from exscope.core.chromosome_coverage import (
        ChromosomeAnalysisResult, 
        ChromosomeCoverage, 
        ChromosomeStatus,
        Karyotype
    )
    
    # Create synthetic result for normal female (XX)
    chromosomes = {}
    autosomes = [f"chr{i}" for i in range(1, 23)]
    
    # Autosomes with ratio ~1.0
    for chrom in autosomes:
        chromosomes[chrom] = ChromosomeCoverage(
            chrom=chrom,
            mean_coverage=150.0,
            median_coverage=145.0,
            percent_covered=98.5,
            total_bases=100000000,
            covered_bases=98500000,
            status=ChromosomeStatus.NORMAL,
            normalized_ratio=1.0
        )
    
    # chrX with ratio ~1.0 (two copies in female)
    chromosomes["chrX"] = ChromosomeCoverage(
        chrom="chrX",
        mean_coverage=150.0,
        median_coverage=145.0,
        percent_covered=98.5,
        total_bases=50000000,
        covered_bases=49250000,
        status=ChromosomeStatus.NORMAL,
        normalized_ratio=1.0
    )
    
    # chrY with ratio ~0.0 (absent in female)
    chromosomes["chrY"] = ChromosomeCoverage(
        chrom="chrY",
        mean_coverage=0.0,
        median_coverage=0.0,
        percent_covered=0.0,
        total_bases=25000000,
        covered_bases=0,
        status=ChromosomeStatus.NOT_SEQUENCED,
        normalized_ratio=0.0
    )
    
    result = ChromosomeAnalysisResult(
        sample_id="test_sample_XX",
        chromosomes=chromosomes,
        autosomal_mean=150.0,
        inferred_karyotype=Karyotype.XX,
        deletions=[],
        warnings=[]
    )
    
    assert result.sample_id == "test_sample_XX"
    assert result.autosomal_mean == 150.0
    assert result.inferred_karyotype == Karyotype.XX
    assert len(result.deletions) == 0
    assert result.chromosomes["chr1"].normalized_ratio == 1.0
    assert result.chromosomes["chrX"].normalized_ratio == 1.0
    assert result.chromosomes["chrY"].normalized_ratio == 0.0


def test_chromosome_analysis_male_karyotype():
    """Test ChromosomeAnalysisResult for normal male (XY)."""
    from exscope.core.chromosome_coverage import (
        ChromosomeAnalysisResult, 
        ChromosomeCoverage, 
        ChromosomeStatus,
        Karyotype
    )
    
    chromosomes = {}
    autosomes = [f"chr{i}" for i in range(1, 23)]
    
    # Autosomes with ratio ~1.0
    for chrom in autosomes:
        chromosomes[chrom] = ChromosomeCoverage(
            chrom=chrom,
            mean_coverage=150.0,
            median_coverage=145.0,
            percent_covered=98.5,
            total_bases=100000000,
            covered_bases=98500000,
            status=ChromosomeStatus.NORMAL,
            normalized_ratio=1.0
        )
    
    # chrX with ratio ~0.5 (one copy in male)
    chromosomes["chrX"] = ChromosomeCoverage(
        chrom="chrX",
        mean_coverage=75.0,
        median_coverage=72.0,
        percent_covered=98.0,
        total_bases=50000000,
        covered_bases=49000000,
        status=ChromosomeStatus.REDUCED,
        normalized_ratio=0.5
    )
    
    # chrY with ratio ~0.5 (one copy in male)
    chromosomes["chrY"] = ChromosomeCoverage(
        chrom="chrY",
        mean_coverage=75.0,
        median_coverage=72.0,
        percent_covered=95.0,
        total_bases=25000000,
        covered_bases=23750000,
        status=ChromosomeStatus.REDUCED,
        normalized_ratio=0.5
    )
    
    result = ChromosomeAnalysisResult(
        sample_id="test_sample_XY",
        chromosomes=chromosomes,
        autosomal_mean=150.0,
        inferred_karyotype=Karyotype.XY,
        deletions=[],
        warnings=[]
    )
    
    assert result.inferred_karyotype == Karyotype.XY
    assert result.chromosomes["chrX"].normalized_ratio == 0.5
    assert result.chromosomes["chrY"].normalized_ratio == 0.5


def test_chromosome_analysis_trisomy21():
    """Test ChromosomeAnalysisResult for trisomy 21 (Down syndrome)."""
    from exscope.core.chromosome_coverage import (
        ChromosomeAnalysisResult, 
        ChromosomeCoverage, 
        ChromosomeStatus,
        Karyotype
    )
    
    chromosomes = {}
    autosomes = [f"chr{i}" for i in range(1, 23)]
    
    # Most autosomes with ratio ~1.0
    for chrom in autosomes:
        if chrom == "chr21":
            # chr21 with ratio ~1.5 (three copies)
            chromosomes[chrom] = ChromosomeCoverage(
                chrom=chrom,
                mean_coverage=225.0,  # 1.5x autosomal mean
                median_coverage=220.0,
                percent_covered=99.0,
                total_bases=50000000,
                covered_bases=49500000,
                status=ChromosomeStatus.ELEVATED,
                normalized_ratio=1.5
            )
        else:
            chromosomes[chrom] = ChromosomeCoverage(
                chrom=chrom,
                mean_coverage=150.0,
                median_coverage=145.0,
                percent_covered=98.5,
                total_bases=100000000,
                covered_bases=98500000,
                status=ChromosomeStatus.NORMAL,
                normalized_ratio=1.0
            )
    
    # Normal female sex chromosomes
    chromosomes["chrX"] = ChromosomeCoverage(
        chrom="chrX",
        mean_coverage=150.0,
        median_coverage=145.0,
        percent_covered=98.5,
        total_bases=50000000,
        covered_bases=49250000,
        status=ChromosomeStatus.NORMAL,
        normalized_ratio=1.0
    )
    
    chromosomes["chrY"] = ChromosomeCoverage(
        chrom="chrY",
        mean_coverage=0.0,
        median_coverage=0.0,
        percent_covered=0.0,
        total_bases=25000000,
        covered_bases=0,
        status=ChromosomeStatus.NOT_SEQUENCED,
        normalized_ratio=0.0
    )
    
    result = ChromosomeAnalysisResult(
        sample_id="test_trisomy21",
        chromosomes=chromosomes,
        autosomal_mean=150.0,
        inferred_karyotype=Karyotype.XX,
        deletions=[],
        warnings=[]
    )
    
    assert result.chromosomes["chr21"].normalized_ratio == 1.5
    assert result.chromosomes["chr21"].status == ChromosomeStatus.ELEVATED


def test_normalized_ratio_plot_creation(tmp_path):
    """Test that plot_normalized_depth_ratio creates valid output."""
    from exscope.visualizers.chromosome_plot import ChromosomePlotter
    from exscope.core.chromosome_coverage import (
        ChromosomeAnalysisResult, 
        ChromosomeCoverage, 
        ChromosomeStatus,
        Karyotype
    )
    
    # Create synthetic result
    chromosomes = {}
    for i in range(1, 23):
        chromosomes[f"chr{i}"] = ChromosomeCoverage(
            chrom=f"chr{i}",
            mean_coverage=150.0,
            median_coverage=145.0,
            percent_covered=98.5,
            total_bases=100000000,
            covered_bases=98500000,
            status=ChromosomeStatus.NORMAL,
            normalized_ratio=1.0
        )
    
    chromosomes["chrX"] = ChromosomeCoverage(
        chrom="chrX",
        mean_coverage=150.0,
        median_coverage=145.0,
        percent_covered=98.5,
        total_bases=50000000,
        covered_bases=49250000,
        status=ChromosomeStatus.NORMAL,
        normalized_ratio=1.0
    )
    
    chromosomes["chrY"] = ChromosomeCoverage(
        chrom="chrY",
        mean_coverage=0.0,
        median_coverage=0.0,
        percent_covered=0.0,
        total_bases=25000000,
        covered_bases=0,
        status=ChromosomeStatus.NOT_SEQUENCED,
        normalized_ratio=0.0
    )
    
    result = ChromosomeAnalysisResult(
        sample_id="test_plot",
        chromosomes=chromosomes,
        autosomal_mean=150.0,
        inferred_karyotype=Karyotype.XX,
        deletions=[],
        warnings=[]
    )
    
    # Create plot
    plotter = ChromosomePlotter()
    output_file = tmp_path / "test_normalized_ratio.html"
    
    try:
        plotter.plot_normalized_depth_ratio(result, output=str(output_file))
        assert output_file.exists()
        
        # Check file size (should have content)
        assert output_file.stat().st_size > 1000
        
        # Check that it's valid HTML
        with open(output_file) as f:
            content = f.read()
            assert "<html>" in content or "<!DOCTYPE html>" in content
            assert "Plotly" in content  # Should use Plotly
            
    except Exception as e:
        pytest.fail(f"Failed to create normalized ratio plot: {e}")


def test_clinical_report_plot_creation(tmp_path):
    """Test that plot_clinical_report creates two-panel output."""
    from exscope.visualizers.chromosome_plot import ChromosomePlotter
    from exscope.core.chromosome_coverage import (
        ChromosomeAnalysisResult, 
        ChromosomeCoverage, 
        ChromosomeStatus,
        Karyotype
    )
    
    # Create synthetic result with deletion
    chromosomes = {}
    for i in range(1, 23):
        if i == 13:
            # chr13 deletion
            chromosomes[f"chr{i}"] = ChromosomeCoverage(
                chrom=f"chr{i}",
                mean_coverage=75.0,
                median_coverage=72.0,
                percent_covered=50.0,
                total_bases=100000000,
                covered_bases=50000000,
                status=ChromosomeStatus.REDUCED,
                normalized_ratio=0.5
            )
        else:
            chromosomes[f"chr{i}"] = ChromosomeCoverage(
                chrom=f"chr{i}",
                mean_coverage=150.0,
                median_coverage=145.0,
                percent_covered=98.5,
                total_bases=100000000,
                covered_bases=98500000,
                status=ChromosomeStatus.NORMAL,
                normalized_ratio=1.0
            )
    
    chromosomes["chrX"] = ChromosomeCoverage(
        chrom="chrX",
        mean_coverage=75.0,
        median_coverage=72.0,
        percent_covered=98.0,
        total_bases=50000000,
        covered_bases=49000000,
        status=ChromosomeStatus.REDUCED,
        normalized_ratio=0.5
    )
    
    chromosomes["chrY"] = ChromosomeCoverage(
        chrom="chrY",
        mean_coverage=75.0,
        median_coverage=72.0,
        percent_covered=95.0,
        total_bases=25000000,
        covered_bases=23750000,
        status=ChromosomeStatus.REDUCED,
        normalized_ratio=0.5
    )
    
    result = ChromosomeAnalysisResult(
        sample_id="test_clinical_report",
        chromosomes=chromosomes,
        autosomal_mean=150.0,
        inferred_karyotype=Karyotype.XY,
        deletions=["chr13"],
        warnings=[]
    )
    
    # Create clinical report plot
    plotter = ChromosomePlotter()
    output_file = tmp_path / "test_clinical_report.html"
    
    try:
        plotter.plot_clinical_report(result, output=str(output_file))
        assert output_file.exists()
        
        # Check file size (should be larger than single panel)
        assert output_file.stat().st_size > 2000
        
        # Check content
        with open(output_file) as f:
            content = f.read()
            assert "<html>" in content or "<!DOCTYPE html>" in content
            assert "Plotly" in content
            # Verify two panels present by checking for two y-axes
            assert "y2" in content  # Second y-axis indicates two panels
            
    except Exception as e:
        pytest.fail(f"Failed to create clinical report plot: {e}")


def test_coverage_breadth_qc_plot(tmp_path):
    """Test that plot_coverage_breadth_qc creates QC metric visualization."""
    from exscope.visualizers.chromosome_plot import ChromosomePlotter
    from exscope.core.chromosome_coverage import (
        ChromosomeAnalysisResult, 
        ChromosomeCoverage, 
        ChromosomeStatus,
        Karyotype
    )
    
    # Create synthetic result
    chromosomes = {}
    for i in range(1, 23):
        chromosomes[f"chr{i}"] = ChromosomeCoverage(
            chrom=f"chr{i}",
            mean_coverage=150.0,
            median_coverage=145.0,
            percent_covered=2.5,  # WES coverage
            total_bases=100000000,
            covered_bases=2500000,
            status=ChromosomeStatus.NORMAL,
            normalized_ratio=1.0
        )
    
    chromosomes["chrX"] = ChromosomeCoverage(
        chrom="chrX",
        mean_coverage=150.0,
        median_coverage=145.0,
        percent_covered=2.8,
        total_bases=50000000,
        covered_bases=1400000,
        status=ChromosomeStatus.NORMAL,
        normalized_ratio=1.0
    )
    
    chromosomes["chrY"] = ChromosomeCoverage(
        chrom="chrY",
        mean_coverage=0.0,
        median_coverage=0.0,
        percent_covered=0.0,
        total_bases=25000000,
        covered_bases=0,
        status=ChromosomeStatus.NOT_SEQUENCED,
        normalized_ratio=0.0
    )
    
    result = ChromosomeAnalysisResult(
        sample_id="test_qc",
        chromosomes=chromosomes,
        autosomal_mean=150.0,
        inferred_karyotype=Karyotype.XX,
        deletions=[],
        warnings=[]
    )
    
    # Create QC plot
    plotter = ChromosomePlotter()
    output_file = tmp_path / "test_qc_breadth.html"
    
    try:
        plotter.plot_coverage_breadth_qc(result, output=str(output_file))
        assert output_file.exists()
        
        # Check content
        with open(output_file) as f:
            content = f.read()
            assert "QC" in content or "Quality" in content
            
    except Exception as e:
        pytest.fail(f"Failed to create QC breadth plot: {e}")
