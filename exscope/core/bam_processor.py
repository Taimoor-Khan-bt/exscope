"""BAM file processing using samtools.

This module provides a wrapper around samtools for efficient BAM file processing.
It uses either subprocess calls to samtools or pysam's built-in functionality.
"""

import subprocess
from pathlib import Path
import logging
import pysam
from typing import Optional, List, Dict, Union
import json

class BamProcessor:
    """Wrapper for samtools operations on BAM files."""
    
    def __init__(self, bam_path: Union[str, Path], reference: Optional[Union[str, Path]] = None, threads: int = 4):
        """Initialize BamProcessor.
        
        Args:
            bam_path: Path to the BAM/CRAM file
            reference: Optional path to reference genome (required for CRAM)
            threads: Number of threads to use for operations
        """
        self.bam_path = Path(bam_path)
        self.reference = Path(reference) if reference else None
        self.threads = threads
        self._validate_files()
        self._check_samtools()
    
    def _validate_files(self):
        """Validate input files existence and format."""
        if not self.bam_path.exists():
            raise FileNotFoundError(f"BAM/CRAM file not found: {self.bam_path}")
        if self.bam_path.suffix.lower() == '.cram' and not self.reference:
            raise ValueError("Reference genome is required for CRAM files")
        if self.reference and not self.reference.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference}")
    
    def _check_samtools(self):
        """Check if samtools is available and get version."""
        try:
            result = subprocess.run(['samtools', '--version'], 
                                 capture_output=True, text=True, check=True)
            logging.info(f"Using samtools: {result.stdout.split(chr(10))[0]}")
        except subprocess.CalledProcessError:
            raise RuntimeError("samtools not found. Please ensure it's installed and in PATH")
    
    def _run_samtools(self, cmd: List[str], capture_output: bool = True) -> subprocess.CompletedProcess:
        """Run a samtools command with proper error handling.
        
        Args:
            cmd: Command list to execute
            capture_output: Whether to capture command output
        """
        try:
            return subprocess.run(cmd, capture_output=capture_output, text=True, check=True)
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if e.stderr else str(e)
            raise RuntimeError(f"Samtools command failed: {error_msg}")
    
    def index(self):
        """Create index for BAM/CRAM file using samtools."""
        cmd = ['samtools', 'index', '-@', str(self.threads), str(self.bam_path)]
        self._run_samtools(cmd)
        
    def sort(self, output: Optional[Path] = None, by_name: bool = False):
        """Sort BAM/CRAM file.
        
        Args:
            output: Output file path (default: input.sorted.bam)
            by_name: Sort by read name instead of coordinates
        """
        if not output:
            output = self.bam_path.with_suffix('.sorted.bam')
            
        cmd = ['samtools', 'sort', '-@', str(self.threads)]
        if by_name:
            cmd.append('-n')
        if self.reference:
            cmd.extend(['--reference', str(self.reference)])
        cmd.extend(['-o', str(output), str(self.bam_path)])
        
        self._run_samtools(cmd)
        return output
    
    def view_region(self, region: str, output: Optional[Path] = None, 
                   min_quality: int = 0, include_flags: Optional[str] = None,
                   exclude_flags: Optional[str] = None) -> Optional[str]:
        """Extract reads from a specific region with filtering options.
        
        Args:
            region: Region in format 'chr:start-end'
            output: Optional output file path
            min_quality: Minimum mapping quality
            include_flags: Required flags (e.g., '0x2' for properly paired)
            exclude_flags: Excluded flags (e.g., '0x4' for unmapped)
        """
        cmd = ['samtools', 'view', '-@', str(self.threads)]
        
        if output:
            cmd.extend(['-b', '-o', str(output)])
        if min_quality > 0:
            cmd.extend(['-q', str(min_quality)])
        if include_flags:
            cmd.extend(['-f', include_flags])
        if exclude_flags:
            cmd.extend(['-F', exclude_flags])
        if self.reference:
            cmd.extend(['--reference', str(self.reference)])
            
        cmd.extend([str(self.bam_path), region])
        
        result = self._run_samtools(cmd)
        if not output:
            return result.stdout
            
    def get_stats(self, output: Optional[Path] = None) -> Dict:
        """Get detailed statistics about the BAM/CRAM file.
        
        Args:
            output: Optional output file path for full stats
        """
        cmd = ['samtools', 'stats', '-@', str(self.threads)]
        if self.reference:
            cmd.extend(['--reference', str(self.reference)])
        cmd.append(str(self.bam_path))
        
        result = self._run_samtools(cmd)
        
        if output:
            output.write_text(result.stdout)
            
        # Parse relevant statistics
        stats = {}
        for line in result.stdout.split('\\n'):
            if line.startswith('SN'):
                fields = line.strip().split('\\t')
                if len(fields) >= 3:
                    key = fields[1].strip(':')
                    value = fields[2]
                    stats[key] = value
                    
        return stats