"""Genomic visualization using pyGenomeTracks.

This module provides interactive visualization capabilities for genomic data
using pyGenomeTracks, with support for multiple track types and customization.
"""

import os
import shutil
import tempfile
from pathlib import Path
import logging
from typing import List, Dict, Optional, Union
from configparser import ConfigParser
from dataclasses import dataclass
import yaml
import matplotlib.pyplot as plt
from pygenometracks import plotTracks

@dataclass
class GenomicRegion:
    """Represents a genomic region for visualization."""
    chrom: str
    start: int
    end: int
    
    def __str__(self) -> str:
        """Convert to string format for pyGenomeTracks."""
        return f"{self.chrom}:{self.start}-{self.end}"

class TrackConfig:
    """Configuration for genome browser tracks."""
    
    TRACK_TYPES = {
        'bigwig': {
            'file_types': ['.bw', '.bigwig'],
            'default_height': 5,
            'default_color': '#33A02C'
        },
        'bed': {
            'file_types': ['.bed', '.bed.gz'],
            'default_height': 3,
            'default_color': '#1F78B4'
        },
        'bedgraph': {
            'file_types': ['.bedgraph', '.bg'],
            'default_height': 4,
            'default_color': '#FF7F00'
        },
        'gtf': {
            'file_types': ['.gtf', '.gff', '.gff3'],
            'default_height': 3,
            'preferred_name': 'gene_name'
        }
    }
    
    def __init__(self, title: str, file_path: Union[str, Path], track_type: str,
                 height: Optional[float] = None, color: Optional[str] = None,
                 **kwargs):
        """Initialize track configuration.
        
        Args:
            title: Track title
            file_path: Path to track data file
            track_type: Type of track (bigwig, bed, bam, gtf)
            height: Track height in the visualization
            color: Track color (for applicable track types)
            **kwargs: Additional track-specific parameters
        """
        self.title = title
        self.file_path = Path(file_path)
        
        if track_type not in self.TRACK_TYPES:
            raise ValueError(f"Unsupported track type: {track_type}")
        self.track_type = track_type
        
        # Validate file format
        if not any(self.file_path.name.endswith(ext) 
                  for ext in self.TRACK_TYPES[track_type]['file_types']):
            raise ValueError(f"Invalid file format for {track_type} track")
            
        # Set default values
        self.height = height or self.TRACK_TYPES[track_type]['default_height']
        self.color = color or self.TRACK_TYPES[track_type].get('default_color')
        self.kwargs = kwargs
        
    def to_dict(self) -> Dict:
        """Convert track configuration to pyGenomeTracks format."""
        config = {
            'file': str(self.file_path),
            'title': self.title,
            'height': self.height
        }
        
        # Add track-specific settings
        if self.track_type == 'bigwig':
            config.update({
                'file_type': 'bigwig',
                'color': self.color,
                'type': 'line'
            })
        elif self.track_type == 'bed':
            config.update({
                'file_type': 'bed',
                'color': self.color,
                'display': 'stacked'
            })
        elif self.track_type == 'bedgraph':
            config.update({
                'file_type': 'bedgraph',
                'color': self.color,
                'type': 'line'
            })
        elif self.track_type == 'gtf':
            config.update({
                'file_type': 'gtf',
                'preferred_name': self.TRACK_TYPES['gtf']['preferred_name'],
                'display': 'stacked'
            })
            
        # Add any additional parameters
        config.update(self.kwargs)
        return config

class GenomeVisualizer:
    """Interactive genome browser using pyGenomeTracks."""
    
    def __init__(self, dpi: int = 100):
        """Initialize visualizer.
        
        Args:
            dpi: Resolution for output images
        """
        self.dpi = dpi
        self.tracks: List[TrackConfig] = []
        # Create a temporary directory for visualization files
        self.config_dir = Path(tempfile.mkdtemp())
        
    def add_track(self, track: TrackConfig):
        """Add a track to the visualization.
        
        Args:
            track: Track configuration
        """
        self.tracks.append(track)
        
    def add_coverage_track(self, file_path: Union[str, Path], title: str = "Coverage",
                          color: str = "#33A02C", height: float = 4):
        """Add a coverage track from BED or bigwig file.
        
        Args:
            file_path: Path to coverage file (BED or bigwig)
            title: Track title
            color: Track color
            height: Track height
        """
        # Detect file format from extension
        path = Path(file_path)
        if path.suffix in ['.bw', '.bigwig']:
            track_type = 'bigwig'
        elif path.name.endswith(('.bed', '.bed.gz')):
            track_type = 'bed'
        else:
            raise ValueError("Coverage file must be in BED or bigwig format")
            
        track = TrackConfig(
            title=title,
            file_path=path,
            track_type=track_type,
            height=height,
            color=color
        )
        self.add_track(track)
        
    def add_annotation_track(self, file_path: Union[str, Path], title: str = "Genes",
                           height: float = 3):
        """Add gene annotation track from GTF/GFF file.
        
        Args:
            file_path: Path to GTF/GFF file
            title: Track title
            height: Track height
        """
        track = TrackConfig(
            title=title,
            file_path=file_path,
            track_type='gtf',
            height=height
        )
        self.add_track(track)
        

        
    def save_config(self, output_path: Union[str, Path]):
        """Save track configuration to file.
        
        Args:
            output_path: Path to save configuration
        """
        config = {}
        for i, track in enumerate(self.tracks, 1):
            config[f'[track{i}]'] = track.to_dict()
            
        with open(output_path, 'w') as f:
            yaml.dump(config, f)
            
    def plot(self, region: GenomicRegion, output: Optional[Union[str, Path]] = None,
            width: float = 40, height_per_track: float = 2):
        """Create visualization for specified region.
        
        Args:
            region: Genomic region to visualize
            output: Optional output file path
            width: Figure width
            height_per_track: Height per track
        """
        if not self.tracks:
            raise ValueError("No tracks added for visualization")
            
        # Create temporary track configuration in config_dir
        config_file = self.config_dir / "tracks.ini"
        config = {}
        
        # Create sections with proper file paths
        for i, track in enumerate(self.tracks, 1):
            section = f'[track{i}]'
            track_config = track.to_dict()
            # Ensure file paths are absolute and properly formatted
            file_path = Path(track_config['file']).resolve()
            track_config['file'] = str(file_path)
            config[section] = track_config
        
        # Write INI configuration
        ini = ConfigParser()
        for section, options in config.items():
            ini.add_section(section.strip('[]'))  # Remove brackets for ConfigParser
            for key, value in options.items():
                ini.set(section.strip('[]'), str(key), str(value))
        
        with open(config_file, 'w') as f:
            ini.write(f)
        
        try:
            # Plot tracks using pygenometracks
            import tempfile
            
            # Get region parameters
            chrom = str(region.chrom)
            start = int(region.start)
            end = int(region.end)
            
            # Use plotTracks function directly
            if output:
                # If output path provided, use it directly
                plotTracks.main([
                    '--tracks', str(config_file),
                    '--region', f'{chrom}:{start}-{end}',
                    '--outFileName', str(output),
                    '--dpi', str(self.dpi)
                ])
                logging.info(f"Visualization saved to: {output}")
            else:
                # Create a temporary file for displaying
                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                    plotTracks.main([
                        '--tracks', str(config_file),
                        '--region', f'{chrom}:{start}-{end}',
                        '--outFileName', tmp.name,
                        '--dpi', str(self.dpi)
                    ])
                    # Display the plot
                    img = plt.imread(tmp.name)
                    plt.imshow(img)
                    plt.axis('off')
                    plt.show()
                    # Clean up
                    os.unlink(tmp.name)
                
        except Exception as e:
            logging.error(f"Failed to create visualization: {str(e)}")
            raise
            
        finally:
            # Clean up temporary directory
            shutil.rmtree(self.config_dir)
            
    def plot_chromosome(self, chromosome: str, output: Union[str, Path],
                       window_size: int = 500000, width: float = 50, height_per_track: float = 2):
        """Create visualization for an entire chromosome.
        
        Args:
            chromosome: Chromosome name (e.g., 'chr17')
            output: Output file path
            window_size: Window size for binning coverage data (default: 500kb)
            width: Figure width
            height_per_track: Height per track
        """
        # For whole chromosome visualization, we need to determine the chromosome length
        # This will typically be obtained from the BAM file header or reference
        # For now, we'll use a large region to capture the entire chromosome
        logging.info(f"Creating visualization for whole {chromosome}")
        
        # Typical chromosome lengths for hg38 (can be made configurable)
        CHROM_LENGTHS = {
            'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
            'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
            'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
            'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
            'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
            'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
            'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
            'chrM': 16569
        }
        
        if chromosome not in CHROM_LENGTHS:
            raise ValueError(f"Unknown chromosome: {chromosome}")
            
        region = GenomicRegion(
            chrom=chromosome,
            start=0,
            end=CHROM_LENGTHS[chromosome]
        )
        
        # Use plot method with adjusted parameters for whole chromosome
        self.plot(region=region, output=output, width=width, height_per_track=height_per_track)
        logging.info(f"Chromosome {chromosome} visualization saved to {output}")
    
    def create_genome_browser(self, region: GenomicRegion,
                          output_html: Union[str, Path]):
        """Create an interactive genome browser visualization.
        
        Args:
            region: Initial genomic region to display
            output_html: Path to save the HTML file
        """
        # This will be implemented when integrating with IGV-reports
        raise NotImplementedError("Interactive genome browser not yet implemented")