"""Tests for BAM processing functionality."""

import pytest
import subprocess
from pathlib import Path
from exscope.core.bam_processor import BamProcessor

def test_bam_processor_initialization(tmp_path):
    """Test BamProcessor initialization."""
    # Should raise error for non-existent file
    with pytest.raises(FileNotFoundError):
        BamProcessor(tmp_path / "nonexistent.bam")
        
    # Create a dummy BAM file
    bam_path = tmp_path / "test.bam"
    bam_path.touch()
    
    # Should initialize successfully
    processor = BamProcessor(bam_path)
    assert processor.bam_path == bam_path
    assert processor.reference is None
    assert processor.threads == 4

def test_bam_processor_cram_validation(tmp_path):
    """Test CRAM file validation."""
    cram_path = tmp_path / "test.cram"
    cram_path.touch()
    
    # Should raise error for CRAM without reference
    with pytest.raises(ValueError, match="Reference genome is required for CRAM files"):
        BamProcessor(cram_path)
        
    # Create dummy reference
    ref_path = tmp_path / "reference.fa"
    ref_path.touch()
    
    # Should initialize successfully with reference
    processor = BamProcessor(cram_path, reference=ref_path)
    assert processor.reference == ref_path

def test_bam_processor_samtools_check(mocker):
    """Test samtools availability check."""
    # Mock successful samtools check
    mocker.patch('subprocess.run', return_value=mocker.Mock(
        stdout="samtools 1.18\\n",
        stderr="",
        returncode=0
    ))
    
    # Should initialize successfully
    bam_path = Path("dummy.bam")
    processor = BamProcessor(bam_path)
    
    # Mock failed samtools check
    mocker.patch('subprocess.run', side_effect=subprocess.CalledProcessError(1, 'cmd'))
    
    # Should raise error
    with pytest.raises(RuntimeError, match="samtools not found"):
        BamProcessor(bam_path)