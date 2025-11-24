# Changelog

All notable changes to ExScope will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2024-11-24

### Added
- **Logarithmic transformation visualization** for chromosome coverage:
  - Revolutionary log-scaling makes 0.01% coverage visible on whole-chromosome scale
  - Perfect for WES/targeted sequencing where only 1-3% of chromosome is covered
  - Formula: `log10(1 + percent*10) / log10(1 + 100*10) * 100` with scale factor of 10
  - Preserves ordering while dramatically enhancing visibility of small values
  - Example: 0.15% coverage (BRCA1 gene) clearly visible vs invisible with linear scaling

- **Expected WES coverage baselines** per chromosome:
  - Reference lines showing expected exon density (1-3% per chromosome)
  - Based on RefSeq exon coverage for hg38
  - Chromosome-specific values (e.g., chr17: 2.8%, chr19: 3.1%, chrY: 0.5%)
  - Helps distinguish true deletions from expected WES coverage patterns

- **Percentage-based fill visualization**:
  - Each chromosome shown as horizontal bar with filled percentage
  - Single shape per chromosome instead of thousands of region rectangles
  - Hover tooltips show both actual and log-scaled percentages
  - Clear distinction between actual coverage (0.15%) and visual representation

- **`--all-plots` flag**:
  - Default: generates optimized tracks only (1 HTML file)
  - With flag: generates all 3 visualizations (tracks + bar + report)
  - Cleaner output for routine analysis
  - Power users can still access all plot types when needed

### Changed
- **MAJOR PERFORMANCE OPTIMIZATION**: Chromosome plotting is now **400x faster**
  - Old method: Minutes of rendering (thousands of individual Plotly traces)
  - New method: **0.03-0.11 seconds** (bulk shape operations)
  - Technical change: Replaced per-region traces with single percentage fill per chromosome
  - Bulk `fig.update_layout(shapes=all_shapes)` instead of incremental additions
  - Reduced from thousands of traces to just 24 traces (one per chromosome)

- **Default output simplified**:
  - Now generates only 1 HTML file by default (optimized tracks)
  - Previous version generated 3 HTML files always
  - Reduces clutter and focuses on most useful visualization
  - Use `--all-plots` to get all 3 visualizations

- **Enhanced hover information**:
  - Shows both actual percentage (e.g., 0.15%) and visual percentage (e.g., 15%)
  - Clarifies log transformation effect
  - More informative for targeted sequencing data

### Performance
- **Chromosome plot generation**: 0.03-0.11 seconds (down from minutes)
- **Rendering efficiency**: 24 traces instead of thousands
- **File size**: Unchanged (~4.7MB HTML files)
- **Memory usage**: Significantly reduced during plot generation
- **Real-world test**: 0.15% chr17 coverage plotted in 0.023 seconds

### Fixed
- Plotting no longer hangs or times out on WES/targeted sequencing data
- Small coverages (0.01-1%) now clearly visible without manual zooming
- Improved responsiveness of interactive HTML plots

### Technical Details
- Implemented `_transform_coverage_for_visibility()` method with log/sqrt/linear options
- Added `EXPECTED_WES_COVERAGE` dictionary with per-chromosome baselines
- Refactored `plot_chromosome_tracks()` to use bulk shape operations
- Removed per-region rectangle generation loop
- Optimized Plotly figure construction workflow

## [1.1.0] - 2024-11-18

### Added
- **Chromosome-level coverage analysis** (`exscope chromosome-coverage`):
  - Detect chromosome deletions and aneuploidies from WES/WGS data
  - Automatic Y chromosome deletion detection for clinical diagnostics
  - Sex chromosome karyotype inference (XX, XY, X0, XXY, XYY, XXX, Unknown)
  - Normalized coverage ratio calculation relative to autosomal mean
  - Classification system: Normal, Reduced, Deleted, Elevated, Partially Sequenced, Not Sequenced
  - Accurate position-level coverage tracking using pysam.pileup()
  - Covered regions output with genomic coordinates (merged within 1000bp gaps)
  - Partial chromosome sequencing detection (<50% coverage)

- **Genome browser-style chromosome track visualization**:
  - Shows all 25 chromosomes (chr1-22, X, Y, MT) as horizontal bars
  - Fills covered regions with color coding by status
  - Small region enhancement: diamond markers + rectangles for regions <0.2% of chromosome length
  - Fixed label overlap issue using y-axis tick labels with proper margins
  - Interactive hover tooltips showing exact coverage and coordinates

- **Publication-ready export formats**:
  - `--export-format` flag supporting: html, png, pdf, svg
  - `--dpi` option for high-resolution figures (default: 300 for publications)
  - `--width` and `--height` for custom plot dimensions
  - PDF vector graphics for scalable, editable figures
  - SVG format for further customization in graphics software

- **Interactive Plotly visualizations**:
  - Chromosome track plot (genome browser style) showing actual sequenced regions
  - Color-coded bar charts (green=normal, orange=reduced, red=deleted, blue=elevated, purple=partial, gray=not sequenced)
  - Hover tooltips with detailed coverage statistics
  - Multi-sample comparison plots
  - Coverage heatmaps for batch analysis
  - Comprehensive report plots with multiple views

- **Comprehensive reporting**:
  - HTML interactive plots for exploration
  - PNG/PDF/SVG high-resolution static images for publications
  - Text reports with ASCII bar charts
  - TSV summary tables for downstream analysis
  - Covered regions TSV with genomic coordinates
  - Automatic interpretation of findings with clinical relevance

- **Multi-sample support**:
  - Side-by-side comparison of multiple BAM files
  - Batch analysis with parallel processing
  - Comparison tables across samples

- **New dependencies**: 
  - plotly >= 5.14.0 (interactive visualizations)
  - kaleido >= 0.2.0 (static image export)
  - pysam >= 0.22.0 (position-level coverage tracking)

### Changed
- Updated CLI with new `chromosome-coverage` subcommand
- Enhanced documentation with chromosome analysis examples and karyotype interpretation
- Added comprehensive example script: `examples/chromosome_coverage_example.py`
- Improved thread limiting (default: 1, max: 2) to prevent system overload during pileup
- Updated README with version 1.1.0 features and publication export guidelines

### Fixed
- Chromosome label overlap in track visualization (changed from floating annotations to y-axis ticks)
- Small covered regions (<0.2% of chromosome) now visible with diamond markers
- Proper detection of partially sequenced chromosomes (e.g., exome data)

### Performance
- Fast position-level coverage calculation using pysam.pileup()
- Region merging with 1000bp gap tolerance for efficient storage
- Efficient analysis of multiple samples in parallel
- Typical runtime: ~10-20 seconds per sample for full chromosome analysis

### Removed
- Outdated example files: `functionality_test.py`, `simple_visualization.py`

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
