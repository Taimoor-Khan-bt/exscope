# ExScope# ExScope# ExScope



ExScope is a comprehensive genomic visualization tool that leverages established bioinformatics algorithms for efficient coverage analysis and publication-quality genome visualizations.



## FeaturesExScope is a comprehensive genomic visualization tool that leverages established bioinformatics algorithms (mosdepth, samtools, pyGenomeTracks) for efficient coverage analysis and beautiful genome-wide visualizations.ExScope is a genomics visualization and analysis tool for visualizing read counts and coverage from BAM/CRAM files for specific gene transcripts, with a focus on strand-specific analysis.



- **Established Tools Integration**: Uses mosdepth for coverage calculation, samtools for BAM/CRAM operations, and pyGenomeTracks for visualization

- **Multi-Mode Visualization**: Single gene, batch processing, or whole chromosome visualization

- **Multiple File Format Support**: BAM/CRAM files with BED/GTF/GFF annotations## Features## Features

- **Parallel Processing**: Fast batch processing with multiprocessing support

- **User-Friendly CLI**: Simple command-line interface with intuitive subcommands



## Installation- **Established Tools Integration**: Uses mosdepth for coverage, samtools for BAM operations, pyGenomeTracks for visualization- Strand-specific coverage visualization



### Prerequisites- **Multi-Mode Visualization**: Single gene, batch processing, or whole chromosome visualization- Support for both BAM and CRAM formats



ExScope requires the following external bioinformatics tools to be installed:- **Multiple File Format Support**: BAM/CRAM files, BED/GTF/GFF annotations- Exon structure visualization with directional indicators



| Tool | Version | Purpose |- **Parallel Processing**: Fast batch processing with multiprocessing- Efficient GTF/GFF parsing with caching system

|------|---------|---------|

| **mosdepth** | >= 0.3.0 | Fast BAM/CRAM depth calculation |- **User-Friendly CLI**: Simple command-line interface with subcommands- Progress tracking for long-running operations

| **samtools** | >= 1.15 | SAM/BAM/CRAM file manipulation |

| **bedtools** | >= 2.30 | Genomic region operations |



#### Install via Conda (Recommended)## Installation## Installation



```bash

# Create a new conda environment (optional but recommended)

conda create -n exscope python=3.8### PrerequisitesTo install ExScope, clone this repository and run the following command from the project's root directory:



# Activate the environment

conda activate exscope

ExScope requires external tools to be installed:```bash

# Install bioinformatics tools

conda install -c bioconda mosdepth samtools bedtools- **mosdepth** (>= 0.3.0): Fast BAM/CRAM depth calculationpip install -e .

```

- **samtools** (>= 1.15): SAM/BAM/CRAM file manipulation```

#### Install via System Package Manager

- **bedtools** (>= 2.30): Genomic region operations

For Ubuntu/Debian:

```bash### Dependencies

sudo apt-get update

sudo apt-get install samtools bedtoolsInstall via conda (recommended):

```

```bashThe tool automatically installs the following dependencies:

For mosdepth, download from [GitHub releases](https://github.com/brentp/mosdepth/releases).

conda install -c bioconda mosdepth samtools bedtools

### Install ExScope

```#### Core Dependencies

1. Clone the repository:

```bash- pysam==0.22.0

git clone https://github.com/Taimoor-Khan-bt/exscope.git

cd exscope### Install ExScope- matplotlib==3.5.1

```

- pandas==1.5.3

2. Install ExScope and Python dependencies:

```bashClone the repository and install:- numpy==1.24.3

pip install -e .

``````bash- gffutils==0.12



This will automatically install all required Python dependencies:git clone <repository-url>- tqdm==4.67.1

- pysam >= 0.22.0

- pandas >= 1.5.0cd exscope

- numpy >= 1.24.0

- gffutils >= 0.12pip install -e .#### Visualization Dependencies

- pybedtools >= 0.9.0

- pygenometracks >= 3.7```- pygenometracks==3.7

- matplotlib >= 3.5.0

- tqdm >= 4.65.0- pybedtools==0.9.1



3. Verify installation:### Dependencies

```bash

exscope --version## Usage

```

#### Core Dependencies

## Usage

- pysam >= 0.22.0```bash

ExScope provides three main operational modes via subcommands:

- pandas >= 1.5.0exscope input.bam annotation.gtf -t TRANSCRIPT_ID -o output.png

### 1. Single Gene Visualization

- numpy >= 1.24.0```

Visualize coverage and annotations for a single gene:

- gffutils >= 0.12

```bash

exscope single <BAM_FILE> <ANNOTATION_FILE> -g <GENE_NAME> -o <OUTPUT_FILE>- pybedtools >= 0.9.0### Arguments

```

- `input.bam`: Input BAM/CRAM file

**Required Arguments:**

- `BAM_FILE`: Input BAM/CRAM file (must be indexed)#### Visualization Dependencies- `annotation.gtf`: Annotation GTF/GFF file

- `ANNOTATION_FILE`: Annotation file in BED, GTF, or GFF format

- `-g, --gene`: Gene name to visualize- pygenometracks >= 3.7- `-t, --transcript`: Transcript ID to visualize

- `-o, --output`: Output file path (PNG/SVG/PDF)

- matplotlib >= 3.5.0- `-o, --output`: Output file (PNG/SVG/PDF)

**Optional Arguments:**

- `--reference`: Reference genome FASTA (required for CRAM files)- `--dpi`: Output image resolution (default: 300)

- `--threads`: Number of threads for processing (default: 4)

- `--dpi`: Image resolution (default: 100)## Usage- `--reference`: Reference genome FASTA (required for CRAM files)



**Example:**- `-v, --version`: Show program's version number and exit

```bash

exscope single sample.bam targets.bed -g BRCA1 -o brca1_coverage.png --threads 4ExScope provides three main modes of operation:- `-h, --help`: Show help message and exit

```



### 2. Batch Processing

### 1. Single Gene Visualization## Example

Process multiple genes in parallel:



```bash

exscope batch <BAM_FILE> <ANNOTATION_FILE> -g <GENE_LIST> -d <OUTPUT_DIR>Visualize coverage and annotations for a single gene:```bash

```

# Using BAM file

**Required Arguments:**

- `BAM_FILE`: Input BAM/CRAM file (must be indexed)```bashexscope sample.bam genes.gtf -t ENST00000357654.9 -o brca1_coverage.png

- `ANNOTATION_FILE`: Annotation file in BED, GTF, or GFF format

- `-g, --genes`: Comma-separated list of gene namesexscope single sample.bam targets.bed -g BRCA1 -o brca1.png

- `-d, --output-dir`: Output directory for results

```# Using CRAM file

**Alternative Input:**

- `--genes-file`: Path to file containing gene names (one per line)exscope sample.cram genes.gtf -t ENST00000357654.9 -o output.png --reference genome.fa



**Optional Arguments:****Arguments:**```

- `--reference`: Reference genome FASTA (required for CRAM files)

- `--threads`: Number of parallel processes (default: 4)- `bam`: Input BAM/CRAM file

- `--no-parallel`: Disable parallel processing (run sequentially)

- `--dpi`: Image resolution (default: 100)- `annotation`: Annotation file (BED/GTF/GFF)## Project Structure



**Examples:**- `-g, --gene`: Gene name to visualize (required)

```bash

# Comma-separated gene list- `-o, --output`: Output file path (required)```

exscope batch sample.bam targets.bed -g "BRCA1,BRCA2,TP53" -d results/ --threads 8

- `--reference`: Reference genome FASTA (for CRAM files)exscope/

# From file

exscope batch sample.bam targets.bed --genes-file genes.txt -d results/- `--threads`: Number of threads (default: 4)├── algorithms/         # Analysis algorithms

```

- `--dpi`: Image resolution (default: 100)│   ├── coverage/      # Coverage calculation

### 3. Chromosome Visualization

│   └── roh/          # Regions of homozygosity

Visualize coverage across an entire chromosome:

### 2. Batch Processing├── core/             # Core functionality

```bash

exscope chromosome <BAM_FILE> <ANNOTATION_FILE> -c <CHROMOSOME> -o <OUTPUT_FILE>│   ├── annotations.py

```

Process multiple genes in parallel:│   ├── bam_processor.py

**Required Arguments:**

- `BAM_FILE`: Input BAM/CRAM file (must be indexed)│   └── coverage.py

- `ANNOTATION_FILE`: Annotation file in BED, GTF, or GFF format

- `-c, --chromosome`: Chromosome name (e.g., chr17, chr1, chrX)```bash└── visualizers/      # Visualization components

- `-o, --output`: Output file path (PNG/SVG/PDF)

# From comma-separated list```

**Optional Arguments:**

- `--reference`: Reference genome FASTA (required for CRAM files)exscope batch sample.bam targets.bed -g "BRCA1,BRCA2,TP53" -d output_dir/

- `--threads`: Number of threads (default: 4)

- `--dpi`: Image resolution (default: 100)## Changelog



**Example:**# From file (one gene per line)

```bash

exscope chromosome sample.bam targets.bed -c chr17 -o chr17_coverage.pngexscope batch sample.bam targets.bed --genes-file genes.txt -d output_dir/### v0.3.0

```

```* Added support for CRAM format

### Global Options

* Added new visualization dependencies (pygenometracks, pybedtools)

Available for all subcommands:

- `--verbose`: Enable detailed logging output**Arguments:*** Improved strand-specific coverage visualization

- `--version`: Show version information and exit

- `--help`: Show help message and exit- `-g, --genes`: Comma-separated list of gene names* Changed CLI to use positional arguments for input files



## Complete Usage Examples- `--genes-file`: File containing gene names (one per line)* Added progress tracking for long operations



### Example 1: Basic Single Gene Visualization- `-d, --output-dir`: Output directory for visualizations (required)* Updated all core dependencies to specific versions

```bash

exscope single Example_deduped.bam xgen-targets-hg38.bed \- `--threads`: Number of parallel processes (default: 4)

  -g BRCA1 \

  -o brca1_visualization.png \- `--no-parallel`: Disable parallel processing### v0.2.0

  --threads 4 \

  --dpi 300* Added a `--version` flag to the command-line tool

```

### 3. Chromosome Visualization* Fixed output file handling

### Example 2: CRAM File with Reference

```bash* Updated version number to 0.2.0

exscope single sample.cram targets.bed \

  -g TP53 \Visualize an entire chromosome:

  -o tp53_coverage.png \

  --reference hg38.fa \### v0.1.0

  --threads 4

``````bash* Initial release



### Example 3: Batch Processing Multiple Genesexscope chromosome sample.bam targets.bed -c chr17 -o chr17.png

```bash

# Create a genes.txt file```## License

echo -e "BRCA1\nBRCA2\nTP53\nPTEN" > genes.txt



# Run batch processing

exscope batch sample.bam targets.bed \**Arguments:**MIT License - See LICENSE file for details.

  --genes-file genes.txt \

  -d batch_results/ \- `-c, --chromosome`: Chromosome name (e.g., chr17)

  --threads 8 \- `-o, --output`: Output file path (required)

  --dpi 150

```### Common Options



### Example 4: Whole Chromosome ViewAll commands support:

```bash- `--verbose`: Enable verbose logging

exscope chromosome wgs.bam targets.bed \- `--version`: Show version information

  -c chr17 \- `--help`: Show help message

  -o chromosome_17_overview.png \

  --threads 4## Examples

```

### Example 1: Single Gene

## Output Files```bash

exscope single ex-data/Example_deduped.bam xgen-targets.bed \

### Single Gene Mode  -g BRCA1 \

- **Visualization file**: High-quality PNG/SVG/PDF with:  -o brca1_coverage.png \

  - Coverage track showing sequencing depth  --threads 4

  - Target region annotations```

  - Gene coordinates and statistics

- **Coverage files**: Mosdepth output files (`.regions.bed.gz`, `.summary.txt`)### Example 2: Batch Processing

```bash

### Batch Modeexscope batch sample.bam targets.bed \

- **Individual visualizations**: One PNG per gene in output directory  -g "BRCA1,BRCA2,TP53,PTEN" \

- **Summary report**: `batch_summary.txt` with:  -d results/ \

  - Processing status for each gene  --threads 8

  - Mean and median coverage statistics```

  - Number of target regions

  - Error messages for failed genes### Example 3: Whole Chromosome

```bash

### Chromosome Modeexscope chromosome wgs.bam targets.bed \

- **Chromosome-wide visualization**: Large-scale coverage plot  -c chr17 \

- **Coverage files**: Whole-chromosome mosdepth outputs  -o chr17_coverage.png \

  --threads 4

## Project Structure```



```### Example 4: CRAM File

exscope/```bash

├── cli.py                    # Command-line interface entry pointexscope single sample.cram targets.bed \

├── core/                     # Core functionality modules  -g BRCA1 \

│   ├── annotations.py        # GTF/GFF/BED annotation handling  -o output.png \

│   ├── bam_processor.py      # BAM/CRAM operations via samtools  --reference hg38.fa

│   ├── batch_processor.py    # Parallel batch processing logic```

│   └── coverage.py           # Coverage calculation via mosdepth

├── visualizers/              # Visualization components## Project Structure

│   └── genomics_track.py     # pyGenomeTracks integration

├── examples/                 # Example scripts```

├── tests/                    # Test suiteexscope/

├── setup.py                  # Installation configuration├── cli.py              # Command-line interface

└── README.md                 # This file├── core/               # Core functionality

```│   ├── annotations.py  # Annotation file handling (GTF/GFF/BED)

│   ├── bam_processor.py # BAM operations via samtools

## Algorithms and Tools│   ├── batch_processor.py # Batch processing with multiprocessing

│   └── coverage.py     # Coverage calculation via mosdepth

ExScope integrates the following established bioinformatics tools:└── visualizers/        # Visualization components

    └── genomics_track.py # pyGenomeTracks integration

### 1. mosdepth```

**Purpose**: Fast BAM/CRAM coverage calculation

## Output

**Citation**: Pedersen BS, Quinlan AR. (2018) *mosdepth: quick coverage calculation for genomes and exomes*. Bioinformatics, 34(5), 867-868. [DOI: 10.1093/bioinformatics/btx699](https://doi.org/10.1093/bioinformatics/btx699)

ExScope generates:

**Features**:- **Single/Batch Mode**: High-quality PNG/SVG/PDF visualizations with:

- Per-base or per-region depth calculation  - Coverage tracks showing depth across gene regions

- 10-20x faster than samtools depth  - Target region annotations

- Multi-threaded processing  - Mean coverage statistics

- Low memory footprint  

- **Batch Mode Summary**: HTML report with:

**GitHub**: [https://github.com/brentp/mosdepth](https://github.com/brentp/mosdepth)  - Success/failure status for each gene

  - Mean coverage per gene

### 2. samtools  - Execution time and errors

**Purpose**: SAM/BAM/CRAM file manipulation

## Algorithms Used

**Citation**: Li H, Handsaker B, Wysoker A, et al. (2009) *The Sequence Alignment/Map format and SAMtools*. Bioinformatics, 25(16), 2078-2079. [DOI: 10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)

1. **mosdepth** ([https://github.com/brentp/mosdepth](https://github.com/brentp/mosdepth))

**Features**:   - Fast, parallelized depth calculation

- Industry-standard BAM/CRAM operations   - Per-base or per-region coverage

- Read extraction and statistics   - Significantly faster than samtools depth

- Sorting, indexing, and format conversion

2. **samtools** ([http://www.htslib.org/](http://www.htslib.org/))

**Website**: [http://www.htslib.org/](http://www.htslib.org/)   - BAM/CRAM file operations

   - Read extraction and statistics

### 3. pyGenomeTracks   - Industry-standard tool

**Purpose**: Publication-quality genome track visualization

3. **pyGenomeTracks** ([https://github.com/deeptools/pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks))

**Citation**: Lopez-Delisle L, Rabbani L, Wolff J, et al. (2021) *pyGenomeTracks: reproducible plots for multivariate genomic datasets*. Bioinformatics, 37(3), 422-423. [DOI: 10.1093/bioinformatics/btaa692](https://doi.org/10.1093/bioinformatics/btaa692)   - Publication-quality genome track visualizations

   - Multiple track types (coverage, genes, variants)

**Features**:   - Highly customizable

- Multiple track types (coverage, genes, variants, etc.)

- Highly customizable visualizations## Performance

- Publication-ready quality

- **mosdepth**: 10-20x faster than samtools depth

**GitHub**: [https://github.com/deeptools/pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)- **Batch Processing**: Linear scaling with CPU cores

- **Memory Efficient**: Streaming operations, low memory footprint

### 4. bedtools

**Purpose**: Genomic interval operations## Changelog



**Citation**: Quinlan AR, Hall IM. (2010) *BEDTools: a flexible suite of utilities for comparing genomic features*. Bioinformatics, 26(6), 841-842. [DOI: 10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)### v1.0.0 (Current)

* Complete refactor to use established bioinformatics tools

**Website**: [https://bedtools.readthedocs.io/](https://bedtools.readthedocs.io/)* New CLI with single/batch/chromosome subcommands

* Integration of mosdepth for coverage calculation

## Performance* Integration of samtools for BAM operations

* Integration of pyGenomeTracks for visualization

- **Coverage Calculation**: 10-20x faster than samtools depth (via mosdepth)* Parallel batch processing support

- **Batch Processing**: Linear scaling with CPU cores using multiprocessing* Comprehensive test suite

- **Memory Usage**: Low memory footprint with streaming operations* Updated documentation

- **Typical Runtime**: 

  - Single gene: ~3-5 seconds### v0.3.0

  - Batch (10 genes): ~30 seconds with 4 threads* Added support for CRAM format

  - Whole chromosome: ~10-15 seconds* Added new visualization dependencies

* Improved strand-specific coverage visualization

## License* Updated core dependencies



MIT License - See LICENSE file for details.### v0.2.0

* Added version flag

## Citation* Fixed output file handling



If you use ExScope in your research, please cite:### v0.1.0

* Initial release

1. **ExScope tool**:

   ```## License

   Khan T. (2024). ExScope: A genomic visualization tool integrating 

   mosdepth, samtools, and pyGenomeTracks. MIT License - See LICENSE file for details.

   GitHub: https://github.com/Taimoor-Khan-bt/exscope

   ```## Citation



2. **Underlying tools** (please cite all):If you use ExScope in your research, please cite the underlying tools:

   - **mosdepth**: Pedersen BS, Quinlan AR. (2018) Bioinformatics, 34(5), 867-868.- mosdepth: Pedersen BS, Quinlan AR. Bioinformatics. 2018

   - **samtools**: Li H, et al. (2009) Bioinformatics, 25(16), 2078-2079.- samtools: Li H, et al. Bioinformatics. 2009

   - **pyGenomeTracks**: Lopez-Delisle L, et al. (2021) Bioinformatics, 37(3), 422-423.- pyGenomeTracks: Lopez-Delisle L, et al. Bioinformatics. 2021

   - **bedtools**: Quinlan AR, Hall IM. (2010) Bioinformatics, 26(6), 841-842.

## Contact

- **Author**: Taimoor Khan
- **Email**: taimoorkhan@scichores.com
- **GitHub**: [https://github.com/Taimoor-Khan-bt/exscope](https://github.com/Taimoor-Khan-bt/exscope)

## Version

Current version: **1.0.0**
