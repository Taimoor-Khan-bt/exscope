# Changelog

All notable changes to ExScope will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-11-18

### Major Refactor
Complete rewrite of ExScope to integrate established bioinformatics tools and provide a professional, modular architecture.

### Added
- **New CLI with three operational modes**:
  - `exscope single`: Single gene visualization
  - `exscope batch`: Batch processing with parallel execution
  - `exscope chromosome`: Whole chromosome visualization
- **Integration of established tools**:
  - mosdepth v0.3.0+ for fast coverage calculation (10-20x faster)
  - samtools v1.15+ for BAM/CRAM operations
  - pyGenomeTracks v3.7+ for publication-quality visualizations
  - bedtools v2.30+ for genomic region operations
- **Modular architecture**:
  - `core/` directory with coverage, annotations, bam_processor, and batch_processor modules
  - `visualizers/` directory with pyGenomeTracks integration
  - `cli.py` as new entry point
- **Parallel batch processing** using multiprocessing for efficient multi-gene analysis
- **Comprehensive test suite** with 7 test modules covering all functionality
- **Example scripts** demonstrating common workflows
- **Global options** for all commands: `--verbose`, `--version`, `--help`
- **Support for CRAM files** with reference genome specification
- **Conda environment specification** (`environment.yml`)

### Changed
- **Complete documentation rewrite** with:
  - Detailed installation guide with conda and system package manager instructions
  - Comprehensive usage examples for all three modes
  - Proper citations with DOIs for all underlying tools
  - Performance benchmarks
  - Project structure documentation
- **Updated dependencies** to use flexible version requirements (>=)
- **Modified setup.py** with new entry point and updated dependencies

### Removed
- **Old monolithic architecture** (`exscope/main.py`)
- **Empty directories** (`algorithms/`, old `cli/`)
- **Redundant documentation files**

### Performance
- 10-20x faster coverage calculation via mosdepth
- Linear scaling with CPU cores for batch processing
- Low memory footprint with streaming operations
- Typical runtimes:
  - Single gene: ~3-5 seconds
  - Batch (10 genes): ~30 seconds with 4 threads
  - Whole chromosome: ~10-15 seconds

### Fixed
- Type checking issues in batch processor for genomic coordinates
- Improved error handling throughout codebase
- Better validation of input files and parameters

## [0.3.0] - 2024

### Added
- Support for CRAM file format
- New visualization dependencies (pygenometracks, pybedtools)
- Progress tracking for long-running operations

### Changed
- Improved strand-specific coverage visualization
- Updated all core dependencies to specific versions
- Changed CLI to use positional arguments for input files

## [0.2.0] - 2024

### Added
- `--version` flag to show program version

### Fixed
- Output file handling issues

### Changed
- Updated version number to 0.2.0

## [0.1.0] - 2024

### Added
- Initial release of ExScope
- Basic genomic visualization functionality
- Strand-specific coverage analysis
- Exon structure visualization with directional indicators
- GTF/GFF parsing with caching system
- Support for BAM format

[1.0.0]: https://github.com/Taimoor-Khan-bt/exscope/releases/tag/v1.0.0
