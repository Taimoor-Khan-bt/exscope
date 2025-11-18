# ExScope

ExScope is a comprehensive genomic visualization tool that leverages established bioinformatics algorithms for efficient coverage analysis and publication-quality genome visualizations.

## Features

- **Established Tools Integration**: Uses mosdepth for coverage calculation, samtools for BAM/CRAM operations, and pyGenomeTracks for visualization
- **Multi-Mode Visualization**: Single gene, batch processing, or whole chromosome visualization
- **Multiple File Format Support**: BAM/CRAM files with BED/GTF/GFF annotations
- **Parallel Processing**: Fast batch processing with multiprocessing support
- **User-Friendly CLI**: Simple command-line interface with intuitive subcommands
- **Strand-specific Coverage Visualization**: Directional coverage analysis
- **Exon Structure Visualization**: Clear exon annotations with directional indicators
- **Efficient GTF/GFF Parsing**: Caching system for fast annotation processing
- **Progress Tracking**: Real-time progress for long-running operations

## Installation

### Prerequisites

ExScope requires the following external bioinformatics tools:

| Tool | Version | Purpose |
|------|---------|---------|
| **mosdepth** | >= 0.3.0 | Fast BAM/CRAM depth calculation |
| **samtools** | >= 1.15 | SAM/BAM/CRAM file manipulation |
| **bedtools** | >= 2.30 | Genomic region operations |

#### Install via Conda (Recommended)

```bash
# Create a new conda environment
conda create -n exscope python=3.8

# Activate the environment
conda activate exscope

# Install bioinformatics tools
conda install -c bioconda mosdepth samtools bedtools
```

#### Install via System Package Manager

For Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install samtools bedtools
```

For mosdepth, download from [GitHub releases](https://github.com/brentp/mosdepth/releases).

### Install ExScope

1. Clone the repository:
```bash
git clone https://github.com/Taimoor-Khan-bt/exscope.git
cd exscope
```

2. Install ExScope and Python dependencies:
```bash
pip install -e .
```

This automatically installs all required Python dependencies:
- pysam >= 0.22.0
- pandas >= 1.5.0
- numpy >= 1.24.0
- gffutils >= 0.12
- pybedtools >= 0.9.0
- pygenometracks >= 3.7
- matplotlib >= 3.5.0
- tqdm >= 4.65.0

3. Verify installation:
```bash
exscope --version
```

## Usage

ExScope provides three main operational modes:

### 1. Single Gene Visualization

Visualize coverage and annotations for a single gene:

```bash
exscope single <BAM_FILE> <ANNOTATION_FILE> -g <GENE_NAME> -o <OUTPUT_FILE>
```

**Required Arguments:**
- `BAM_FILE`: Input BAM/CRAM file (must be indexed)
- `ANNOTATION_FILE`: Annotation file in BED, GTF, or GFF format
- `-g, --gene`: Gene name to visualize
- `-o, --output`: Output file path (PNG/SVG/PDF)

**Optional Arguments:**
- `--reference`: Reference genome FASTA (required for CRAM files)
- `--threads`: Number of threads (default: 4)
- `--dpi`: Image resolution (default: 100)

**Example:**
```bash
exscope single sample.bam targets.bed -g BRCA1 -o brca1_coverage.png --threads 4
```

### 2. Batch Processing

Process multiple genes in parallel:

```bash
exscope batch <BAM_FILE> <ANNOTATION_FILE> -g <GENE_LIST> -d <OUTPUT_DIR>
```

**Required Arguments:**
- `BAM_FILE`: Input BAM/CRAM file (must be indexed)
- `ANNOTATION_FILE`: Annotation file in BED, GTF, or GFF format
- `-g, --genes`: Comma-separated list of gene names
- `-d, --output-dir`: Output directory for results

**Alternative Input:**
- `--genes-file`: Path to file containing gene names (one per line)

**Optional Arguments:**
- `--reference`: Reference genome FASTA (required for CRAM files)
- `--threads`: Number of parallel processes (default: 4)
- `--no-parallel`: Disable parallel processing
- `--dpi`: Image resolution (default: 100)

**Examples:**
```bash
# Comma-separated gene list
exscope batch sample.bam targets.bed -g "BRCA1,BRCA2,TP53" -d results/ --threads 8

# From file
exscope batch sample.bam targets.bed --genes-file genes.txt -d results/
```

### 3. Chromosome Visualization

Visualize coverage across an entire chromosome:

```bash
exscope chromosome <BAM_FILE> <ANNOTATION_FILE> -c <CHROMOSOME> -o <OUTPUT_FILE>
```

**Required Arguments:**
- `BAM_FILE`: Input BAM/CRAM file (must be indexed)
- `ANNOTATION_FILE`: Annotation file in BED, GTF, or GFF format
- `-c, --chromosome`: Chromosome name (e.g., chr17, chr1, chrX)
- `-o, --output`: Output file path (PNG/SVG/PDF)

**Optional Arguments:**
- `--reference`: Reference genome FASTA (required for CRAM files)
- `--threads`: Number of threads (default: 4)
- `--dpi`: Image resolution (default: 100)

**Example:**
```bash
exscope chromosome sample.bam targets.bed -c chr17 -o chr17_coverage.png
```

### 4. Chromosome-Level Coverage Analysis (NEW)

Analyze coverage across all chromosomes to detect deletions and chromosomal abnormalities:

```bash
exscope chromosome-coverage <BAM_FILE> -d <OUTPUT_DIR>
```

**Required Arguments:**
- `BAM_FILE`: Input BAM/CRAM file (single sample mode)
- `-d, --output-dir`: Output directory for results

**Optional Arguments:**
- `--bams`: Multiple BAM files for comparison (e.g., `--bams sample1.bam sample2.bam`)
- `--bams-file`: File containing BAM paths (one per line)
- `--reference`: Reference genome FASTA (required for CRAM files)
- `--threads`: Number of threads (default: 1, max: 2 to prevent system overload)
- `--dpi`: Image resolution for plots (default: 300 for publication quality)
- `--export-format`: Output format(s): `html`, `png`, `pdf`, `svg` (default: `html`)
- `--width`: Plot width in pixels (default: 1400)
- `--height`: Plot height in pixels (default: 800)

**Output Files:**
- `<sample>_chromosome_tracks.html`: Genome browser-style visualization showing actual sequenced regions
- `<sample>_chromosome_coverage.html`: Interactive bar chart with hover tooltips
- `<sample>_chromosome_report.html`: Comprehensive report with multiple views
- `<sample>_chromosome_report.txt`: Human-readable text report with ASCII bars
- `<sample>_chromosome_summary.tsv`: Tab-delimited coverage statistics
- `<sample>_covered_regions.tsv`: Genomic coordinates of all covered regions
- `chromosome_comparison.html`: Multi-sample comparison (when multiple BAMs provided)
- `chromosome_heatmap.html`: Coverage heatmap across samples

**Single Sample Example:**
```bash
exscope chromosome-coverage patient.bam -d chr_analysis/
```

**Publication-Ready Figures (High-Resolution PDF/PNG):**
```bash
# Generate PNG and PDF at 300 DPI for publications
exscope chromosome-coverage patient.bam -d figures/ \
  --export-format png pdf \
  --dpi 300 \
  --width 1600 \
  --height 900

# Output: *.png and *.pdf publication-ready figures
```

**Multi-Sample Comparison:**
```bash
exscope chromosome-coverage --bams normal.bam patient.bam -d comparison/
```

**Interactive HTML (for exploration):**
```bash
# Default: generates interactive HTML plots
exscope chromosome-coverage patient.bam -d results/ --export-format html
```

**Use Case: Y Chromosome Deletion Detection**

ExScope automatically detects Y chromosome deletions in WES data:

```bash
exscope chromosome-coverage sample.bam -d y_deletion_check/
```

The tool will:
1. Calculate coverage for all chromosomes (chr1-22, X, Y, MT) using pileup analysis
2. Track actual genomic regions with coverage (merged within 1000bp gaps)
3. Normalize coverage to autosomal mean
4. Flag chromosomes with <15% expected coverage as deleted
5. Infer sex chromosome karyotype (XX, XY, X0, XXY, XYY, XXX)
6. Generate genome browser-style visualizations showing covered vs uncovered regions
7. Create interactive plots and detailed reports

**Visualization Features:**
- **Chromosome Track Plot**: Shows all 25 chromosomes as horizontal bars with actual covered regions highlighted
- **Small Region Enhancement**: Regions <0.2% of chromosome length shown with diamond markers for visibility
- **Color Coding**: Normal (green), Reduced (orange), Deleted (red), Elevated (blue), Partial (purple), Not Sequenced (gray)
- **Interactive HTML**: Hover tooltips showing exact coverage values and genomic coordinates

**Karyotype Inference:**
ExScope automatically infers karyotype from sex chromosome coverage ratios:
- **XX** (Normal female): chrX ratio ~1.0, chrY ratio <0.15
- **XY** (Normal male): chrX ratio ~0.5, chrY ratio ~0.5
- **X0** (Turner syndrome): chrX ratio ~0.5, chrY ratio <0.15
- **XXY** (Klinefelter): chrX ratio ~1.0, chrY ratio ~0.5
- **XYY** syndrome: chrX ratio ~0.5, chrY ratio ~1.0
- **XXX** (Triple X): chrX ratio ~1.0, chrY ratio <0.15
- **Unknown**: When ratios don't match expected patterns (e.g., targeted sequencing without X/Y coverage)

**Interpretation:**
- **Normal ratio** (~1.0 for autosomes): Chromosome present with expected copy number
- **Reduced ratio** (0.4-0.6 for X/Y in males): Expected for sex chromosomes  
- **Deleted** (<0.15): Chromosome likely absent or severely underrepresented
- **Elevated** (>1.45): Possible duplication or trisomy
- **Partially Sequenced** (<50% coverage): Targeted sequencing (e.g., exome) with incomplete chromosome coverage

### Global Options

Available for all subcommands:
- `--verbose`: Enable detailed logging
- `--version`: Show version information
- `--help`: Show help message

## Complete Examples

### Example 1: Basic Single Gene
```bash
exscope single Example_deduped.bam xgen-targets-hg38.bed \
  -g BRCA1 \
  -o brca1_visualization.png \
  --threads 4 \
  --dpi 300
```

### Example 2: CRAM File with Reference
```bash
exscope single sample.cram targets.bed \
  -g TP53 \
  -o tp53_coverage.png \
  --reference hg38.fa \
  --threads 4
```

### Example 3: Batch Processing
```bash
# Create a genes.txt file
echo -e "BRCA1\nBRCA2\nTP53\nPTEN" > genes.txt

# Run batch processing
exscope batch sample.bam targets.bed \
  --genes-file genes.txt \
  -d batch_results/ \
  --threads 8
```

### Example 4: Whole Chromosome
```bash
exscope chromosome wgs.bam targets.bed \
  -c chr17 \
  -o chromosome_17.png \
  --threads 4
```

## Output Files

### Single Gene Mode
- **Visualization**: High-quality PNG/SVG/PDF with coverage tracks and annotations
- **Coverage files**: Mosdepth outputs (`.regions.bed.gz`, `.summary.txt`)

### Batch Mode
- **Individual visualizations**: One PNG per gene in output directory
- **Summary report**: `batch_summary.txt` with processing status and statistics

### Chromosome Mode
- **Chromosome-wide visualization**: Large-scale coverage plot
- **Coverage files**: Whole-chromosome mosdepth outputs

## Project Structure

```
exscope/
├── cli.py                    # Command-line interface
├── core/                     # Core functionality
│   ├── annotations.py        # GTF/GFF/BED handling
│   ├── bam_processor.py      # BAM/CRAM operations
│   ├── batch_processor.py    # Batch processing
│   └── coverage.py           # Coverage calculation
├── visualizers/              # Visualization
│   └── genomics_track.py     # pyGenomeTracks integration
├── examples/                 # Example scripts
├── tests/                    # Test suite
└── setup.py                  # Installation config
```

## Algorithms and Tools

ExScope integrates the following established bioinformatics tools:

### 1. mosdepth
**Purpose**: Fast BAM/CRAM coverage calculation  
**Citation**: Pedersen BS, Quinlan AR. (2018) *mosdepth: quick coverage calculation for genomes and exomes*. Bioinformatics, 34(5), 867-868. [DOI: 10.1093/bioinformatics/btx699](https://doi.org/10.1093/bioinformatics/btx699)  
**Features**: Per-base or per-region depth calculation, 10-20x faster than samtools depth, multi-threaded processing, low memory footprint  
**GitHub**: [https://github.com/brentp/mosdepth](https://github.com/brentp/mosdepth)

### 2. samtools
**Purpose**: SAM/BAM/CRAM file manipulation  
**Citation**: Li H, Handsaker B, Wysoker A, et al. (2009) *The Sequence Alignment/Map format and SAMtools*. Bioinformatics, 25(16), 2078-2079. [DOI: 10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)  
**Features**: Industry-standard BAM/CRAM operations, read extraction and statistics, sorting, indexing, and format conversion  
**Website**: [http://www.htslib.org/](http://www.htslib.org/)

### 3. pyGenomeTracks
**Purpose**: Publication-quality genome track visualization  
**Citation**: Lopez-Delisle L, Rabbani L, Wolff J, et al. (2021) *pyGenomeTracks: reproducible plots for multivariate genomic datasets*. Bioinformatics, 37(3), 422-423. [DOI: 10.1093/bioinformatics/btaa692](https://doi.org/10.1093/bioinformatics/btaa692)  
**Features**: Multiple track types (coverage, genes, variants), highly customizable visualizations, publication-ready quality  
**GitHub**: [https://github.com/deeptools/pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)

### 4. bedtools
**Purpose**: Genomic interval operations  
**Citation**: Quinlan AR, Hall IM. (2010) *BEDTools: a flexible suite of utilities for comparing genomic features*. Bioinformatics, 26(6), 841-842. [DOI: 10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)  
**Website**: [https://bedtools.readthedocs.io/](https://bedtools.readthedocs.io/)

## Performance

- **Coverage Calculation**: 10-20x faster than samtools depth (via mosdepth)
- **Batch Processing**: Linear scaling with CPU cores
- **Memory Usage**: Low memory footprint with streaming operations
- **Typical Runtime**:
  - Single gene: ~3-5 seconds
  - Batch (10 genes): ~30 seconds with 4 threads
  - Whole chromosome: ~10-15 seconds

## License

MIT License - See LICENSE file for details.

## Citation

If you use ExScope in your research, please cite:

1. **ExScope tool**:
   ```
   Khan T. (2024). ExScope: A genomic visualization tool integrating 
   mosdepth, samtools, and pyGenomeTracks. 
   GitHub: https://github.com/Taimoor-Khan-bt/exscope
   ```

2. **Underlying tools** (please cite all):
   - **mosdepth**: Pedersen BS, Quinlan AR. (2018) Bioinformatics, 34(5), 867-868.
   - **samtools**: Li H, et al. (2009) Bioinformatics, 25(16), 2078-2079.
   - **pyGenomeTracks**: Lopez-Delisle L, et al. (2021) Bioinformatics, 37(3), 422-423.
   - **bedtools**: Quinlan AR, Hall IM. (2010) Bioinformatics, 26(6), 841-842.

## Contact

- **Author**: Taimoor Khan
- **Email**: taimoorkhan@scichores.com
- **GitHub**: [https://github.com/Taimoor-Khan-bt/exscope](https://github.com/Taimoor-Khan-bt/exscope)

## Version History

- **1.1.0** (2024-11-18):
  - Added chromosome coverage analysis with region tracking
  - Genome browser-style chromosome track visualization
  - Karyotype inference (XX, XY, X0, XXY, XYY, XXX)
  - Publication-ready exports (PNG, PDF, SVG at 300+ DPI)
  - Small region visibility enhancement with diamond markers
  - Covered regions output (genomic coordinates)
  - Partially sequenced chromosome detection

- **1.0.0** (2024-11-17):
  - Initial release
  - Single gene, batch, and chromosome visualization modes
  - Integration with mosdepth, samtools, pyGenomeTracks

Current version: **1.1.0**
