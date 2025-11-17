"""Annotation file processing using gffutils and bedtools.

This module provides comprehensive functionality for handling genomic annotations,
with support for multiple formats, validation, and advanced operations using bedtools.
"""

import gffutils
from pathlib import Path
import subprocess
import logging
import pandas as pd
from typing import Optional, Union, List, Dict, Set
from enum import Enum
import pybedtools
import tempfile
import shutil
import os

class AnnotationFormat(Enum):
    """Supported annotation file formats."""
    GTF = "gtf"
    GFF = "gff"
    GFF3 = "gff3"
    BED = "bed"
    BEDGRAPH = "bedgraph"
    
class GenomicFeature(Enum):
    """Common genomic feature types."""
    GENE = "gene"
    TRANSCRIPT = "transcript"
    EXON = "exon"
    CDS = "CDS"
    UTR = "UTR"
    PROMOTER = "promoter"

class AnnotationProcessor:
    """Handler for genomic annotations with advanced processing capabilities."""
    
    def __init__(self, annotation_path: Union[str, Path]):
        """Initialize AnnotationProcessor.
        
        Args:
            annotation_path: Path to the annotation file (GTF/GFF/BED)
        """
        self.annotation_path = Path(annotation_path)
        self.format = self._detect_format()
        self._validate_file()
        
        # Initialize database if it's a GTF/GFF file
        self.db = None
        # Initialize bedtools
        self.bedtool = None
        
        self._init_bedtools()
        
        if self.format in [AnnotationFormat.GTF, AnnotationFormat.GFF, AnnotationFormat.GFF3]:
            self.create_gff_db()
    
    def _detect_format(self) -> AnnotationFormat:
        """Detect the format of the annotation file."""
        suffix = self.annotation_path.suffix.lower()
        format_map = {
            '.gtf': AnnotationFormat.GTF,
            '.gff': AnnotationFormat.GFF,
            '.gff3': AnnotationFormat.GFF3,
            '.bed': AnnotationFormat.BED,
            '.bedgraph': AnnotationFormat.BEDGRAPH
        }
        if suffix not in format_map:
            raise ValueError(f"Unsupported file format: {suffix}")
        return format_map[suffix]
    
    def _validate_file(self):
        """Validate input file existence and basic format."""
        if not self.annotation_path.exists():
            raise FileNotFoundError(f"Annotation file not found: {self.annotation_path}")
        
        # Check if file is not empty
        if self.annotation_path.stat().st_size == 0:
            raise ValueError(f"Annotation file is empty: {self.annotation_path}")
            
        # Basic format validation
        with open(self.annotation_path) as f:
            first_line = f.readline().strip()
            if self.format in [AnnotationFormat.GTF, AnnotationFormat.GFF, AnnotationFormat.GFF3]:
                if not first_line.startswith('#') and len(first_line.split('\t')) != 9:
                    raise ValueError(f"Invalid {self.format.value} format")
            elif self.format == AnnotationFormat.BED and len(first_line.split('\t')) < 3:
                raise ValueError("Invalid BED format")
    
    def _init_bedtools(self):
        """Initialize pybedtools interface."""
        try:
            if not self.annotation_path.exists():
                raise FileNotFoundError(f"Annotation file not found: {self.annotation_path}")
                
            # Convert to absolute path
            abs_path = self.annotation_path.absolute()
            logging.info(f"Initializing bedtools with file: {abs_path}")
            
            # Create BedTool instance
            self.bedtool = pybedtools.BedTool(str(abs_path))
            if not self.bedtool:
                raise RuntimeError("BedTool initialization failed")
                
            # Verify the file is readable
            _ = list(self.bedtool[0:1])  # Try to read the first line
            logging.info("BedTool initialized successfully")
            
        except Exception as e:
            logging.error(f"Failed to initialize bedtools: {str(e)}")
            raise RuntimeError(f"Failed to initialize bedtools: {str(e)}")
    
    def create_gff_db(self, force: bool = False, merge_strategy: str = "create_unique"):
        """Create or load gffutils database with advanced options.
        
        Args:
            force: Whether to force database recreation
            merge_strategy: Strategy for handling duplicate features
        """
        if self.format not in [AnnotationFormat.GTF, AnnotationFormat.GFF, AnnotationFormat.GFF3]:
            raise ValueError("GFF database can only be created for GTF/GFF files")
            
        db_path = f"{self.annotation_path}.db"
        if force or not Path(db_path).exists():
            self.db = gffutils.create_db(
                str(self.annotation_path),
                dbfn=db_path,
                force=True,
                keep_order=True,
                merge_strategy=merge_strategy,
                disable_infer_genes=False,
                disable_infer_transcripts=False
            )
        else:
            self.db = gffutils.FeatureDB(db_path)
    
    def extract_features(self, feature_type: Union[GenomicFeature, str], 
                        output: Path,
                        upstream: int = 0,
                        downstream: int = 0) -> Path:
        """Extract specific features and optionally include flanking regions.
        
        Args:
            feature_type: Type of feature to extract
            output: Output file path
            upstream: Bases to include upstream
            downstream: Bases to include downstream
        """
        if isinstance(feature_type, GenomicFeature):
            feature_type = feature_type.value
            
        if not self.bedtool:
            raise RuntimeError("BedTool not initialized")
            
        # Create temporary file for filtered features
        with tempfile.NamedTemporaryFile(suffix='.bed') as filtered_tmp:
            if self.format in [AnnotationFormat.GTF, AnnotationFormat.GFF, AnnotationFormat.GFF3]:
                # Ensure GFF database is initialized
                if not self.db:
                    self.create_gff_db()
                if not self.db:
                    raise RuntimeError("Failed to initialize GFF database")
                
                # Convert GFF features to BED format
                with open(filtered_tmp.name, 'w') as f:
                    try:
                        for feature in self.db.features_of_type(feature_type):
                            bed_line = f"{feature.seqid}\t{feature.start}\t{feature.end}\t"
                            bed_line += f"{feature.id}\t{feature.score if feature.score != '.' else '0'}\t{feature.strand}\n"
                            f.write(bed_line)
                    except Exception as e:
                        raise RuntimeError(f"Failed to extract features from GFF database: {str(e)}")
            else:
                # For BED files, filter by column 6 using awk
                filter_cmd = f"awk '$6==\"{feature_type}\"' {self.bedtool.fn}"
                result = subprocess.run(
                    filter_cmd, shell=True, capture_output=True, text=True, check=True)
                with open(filtered_tmp.name, 'w') as f:
                    f.write(result.stdout)
            
            # Create BedTool from filtered features
            filtered = pybedtools.BedTool(filtered_tmp.name)
            
            if upstream > 0 or downstream > 0:
                # Add slop using bedtools command
                with tempfile.NamedTemporaryFile(suffix='.bed') as slop_tmp:
                    with open(slop_tmp.name, 'w') as outfile:
                        result = subprocess.run(
                            [str(p) for p in ['bedtools', 'slop', '-i', filtered.fn, 
                             '-l', str(upstream), '-r', str(downstream)]],
                            capture_output=True, text=True, check=True)
                        outfile.write(result.stdout)
                    filtered = pybedtools.BedTool(slop_tmp.name)
                
            # Save final result
            filtered.saveas(str(output))
            
        return output
    
    def find_overlaps(self, other: Union[str, Path], 
                     min_overlap: float = 0.0,
                     strand_specific: bool = False,
                     invert: bool = False) -> pybedtools.BedTool:
        """Find overlapping regions between two annotation files.
        
        Args:
            other: Path to other annotation file
            min_overlap: Minimum overlap fraction required
            strand_specific: Consider strand when finding overlaps
            invert: Return non-overlapping regions instead
        """
        if not self.bedtool:
            raise RuntimeError("BedTool not initialized")
            
        try:
            # Create BedTool from other file
            other_bed = pybedtools.BedTool(str(other))
            
            # Validate and adjust min_overlap
            if not 0 <= min_overlap <= 1:
                raise ValueError("min_overlap must be between 0 and 1")
            
            # Create output file path
            output_file = self.annotation_path.parent / f"intersect_{os.getpid()}.bed"
            
            # Build command
            cmd = ['bedtools', 'intersect',
                  '-a', str(self.bedtool.fn),
                  '-b', str(other_bed.fn),
                  '-f', str(max(min_overlap, 1e-9))]
            
            if strand_specific:
                cmd.append('-s')
            if invert:
                cmd.append('-v')
                
            # Run command and write output
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            with open(output_file, 'w') as f:
                f.write(result.stdout)
                
            # Return new BedTool
            result = pybedtools.BedTool(str(output_file))
            
            return result
            
        except Exception as e:
            raise RuntimeError(f"Failed to find overlaps: {str(e)}")
    
    def merge_features(self, distance: int = 0, 
                      strand_specific: bool = False) -> pybedtools.BedTool:
        """Merge overlapping or nearby features.
        
        Args:
            distance: Maximum distance between features to merge
            strand_specific: Consider strand when merging
        """
        if not self.bedtool:
            raise RuntimeError("BedTool not initialized")
            
        sorted_bed = self.bedtool.sort()
        if not sorted_bed:
            raise RuntimeError("Failed to sort BedTool")
            
        result = sorted_bed.merge(
            d=distance,
            s=strand_specific,
            c="4,5,6",  # Columns to report for merged features
            o="distinct,mean,distinct"  # Operations to use when merging
        )
        return result
    
    def get_coverage(self, other: Union[str, Path]) -> pd.DataFrame:
        """Calculate coverage statistics for regions.
        
        Args:
            other: Path to BAM/BED file to calculate coverage from
        """
        if not self.bedtool:
            raise RuntimeError("BedTool not initialized")
        
        other_bed = pybedtools.BedTool(str(other))
        
        # Create a temporary file for output
        with tempfile.NamedTemporaryFile(suffix='.bed') as tmp:
            # Run bedtools coverage command
            cmd = [str(p) for p in ['bedtools', 'coverage', '-a', self.bedtool.fn, '-b', other_bed.fn, '-counts']]
            with open(tmp.name, 'w') as outfile:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                outfile.write(result.stdout)
            
            # Read the results
            return pd.read_table(tmp.name, header=None)
    
    def to_format(self, output: Path, format: AnnotationFormat):
        """Convert annotation to another format.
        
        Args:
            output: Output file path
            format: Desired output format
        """
        if not self.bedtool:
            raise RuntimeError("BedTool not initialized")
            
        if format == AnnotationFormat.BED and self.format != AnnotationFormat.BED:
            if self.format in [AnnotationFormat.GTF, AnnotationFormat.GFF, AnnotationFormat.GFF3]:
                # Convert each feature to BED format
                converted = pybedtools.BedTool([
                    self._gtf_to_bed_converter(feature) 
                    for feature in self.bedtool
                ])
                converted.saveas(str(output))
        elif format in [AnnotationFormat.GTF, AnnotationFormat.GFF, AnnotationFormat.GFF3]:
            # Implement conversion to GTF/GFF
            raise NotImplementedError(f"Conversion to {format.value} not yet implemented")
    
    @staticmethod
    def _gtf_to_bed_converter(feature):
        """Convert GTF/GFF feature to BED format."""
        return pybedtools.create_interval_from_list([
            feature.chrom,
            str(feature.start),
            str(feature.end),
            feature.name if feature.name != '.' else feature.fields[2],
            str(feature.score if feature.score != '.' else 0),
            feature.strand
        ])
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        pybedtools.cleanup()  # Clean up any temporary files