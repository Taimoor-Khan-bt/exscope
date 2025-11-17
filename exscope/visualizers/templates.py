"""Example templates for genomic visualization."""

from pathlib import Path
from typing import Optional, Union
from .genomics_track import GenomeVisualizer, TrackConfig, GenomicRegion

def create_coverage_view(coverage_file: Union[str, Path],
                        region: GenomicRegion,
                        output: Optional[Union[str, Path]] = None):
    """Create a basic coverage view.
    
    Args:
        coverage_file: Path to bigwig coverage file
        region: Genomic region to visualize
        output: Optional output file path
    """
    vis = GenomeVisualizer()
    vis.add_coverage_track(coverage_file)
    vis.plot(region, output)

def create_genes_and_coverage(coverage_file: Union[str, Path],
                            annotation_file: Union[str, Path],
                            region: GenomicRegion,
                            output: Optional[Union[str, Path]] = None):
    """Create a view with gene annotations and coverage.
    
    Args:
        coverage_file: Path to bigwig coverage file
        annotation_file: Path to GTF/GFF annotation file
        region: Genomic region to visualize
        output: Optional output file path
    """
    vis = GenomeVisualizer()
    vis.add_coverage_track(coverage_file)
    vis.add_annotation_track(annotation_file)
    vis.plot(region, output)

def create_full_view(coverage_file: Union[str, Path],
                    bam_file: Union[str, Path],
                    annotation_file: Union[str, Path],
                    region: GenomicRegion,
                    output: Optional[Union[str, Path]] = None):
    """Create a comprehensive view with coverage, alignments and annotations.
    
    Args:
        coverage_file: Path to bigwig coverage file
        bam_file: Path to BAM alignment file
        annotation_file: Path to GTF/GFF annotation file
        region: Genomic region to visualize
        output: Optional output file path
    """
    vis = GenomeVisualizer()
    vis.add_coverage_track(coverage_file)
    vis.add_alignment_track(bam_file)
    vis.add_annotation_track(annotation_file)
    vis.plot(region, output)