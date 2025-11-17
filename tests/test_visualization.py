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