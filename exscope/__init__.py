"""ExScope: A comprehensive genomic visualization tool."""

from .core.annotations import AnnotationProcessor, AnnotationFormat
from .core.coverage import CoverageCalculator, SequencingType
from .core.bam_processor import BamProcessor
from .visualizers.genomics_track import GenomeVisualizer, TrackConfig, GenomicRegion

__version__ = "0.3.0"

__all__ = [
    'AnnotationProcessor',
    'AnnotationFormat',
    'CoverageCalculator',
    'SequencingType',
    'BamProcessor',
    'GenomeVisualizer',
    'TrackConfig',
    'GenomicRegion'
]