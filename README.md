# ExScope - Visualizing Read Counts at Genomic Regions

ExScope is a Python-based bioinformatics tool that enables users to visualize read counts at specific genomic regions and Ensembl transcript IDs. It normalizes the read counts and generates plots. This tool is handy for identifying copy number variations (CNVs) and analyzing gene expression of a specific region on chromosome.

## Features
- **Read Count Visualization**: Generates plots showing read counts across specified genomic regions.
- **Normalization**: Normalizes read counts for accurate comparisons.
- **Ensembl Transcript ID Support**: Allows focused analysis on specific Ensembl transcript IDs.

## Installation

To install ExScope, follow these steps:

1. Clone the repository:

    ```bash
    git clone https://github.com/Taimoor-Khan-bt/exscope.git
    cd exscope
    ```

2. Install the package with `pip`:

    ```bash
    pip install .
    ```

This will automatically install all dependencies specified in the `setup.py` file.

## Prerequisites

ExScope requires the following tools and libraries:
- Python 3.6+
- `pysam` (for reading BAM files)
- `matplotlib` (for plotting)
- `pandas` (for data manipulation)
- `scipy` (for clustering and dendrogram generation)
- `argparse` (for command-line argument parsing)

These dependencies will be installed during the installation process.

## Input Files

### 1. BAM File
You will need a BAM file containing aligned sequencing reads as input.

### 2. GFF3 or GTF File
To provide exon annotations, you will need a GFF3 or GTF file. These files can be downloaded from Ensembl, UCSC Genome Browser, or other genomic databases. Ensure the file matches the reference genome used in your analysis.

- **Ensembl GFF3 files**: [Ensembl FTP](https://ftp.ensembl.org/pub/current_gff3/)
- **UCSC GTF files**: [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/)

### 3. Genomic Region and Transcript ID
Specify the genomic region in the format `chr:start-end` and provide the Ensembl transcript ID you wish to analyze.

## Usage

Once installed, you can run ExScope from the command line. Basic usage:

```bash
exscope -b /path/to/your.bam -g /path/to/annotations.gff3 -r chr1:1000000-1050000 -tid ENST00000367770 -o /path/to/output_dir
```

### Command-Line Arguments
- `-b, --bam`: Path to the input BAM file (required).
- `-g, --gff3`: Path to the GFF3 file for exon annotations (required).
- `-r, --region`: Genomic region in the format `chr:start-end` (required).
- `-tid, --transcript_id`: Ensembl transcript ID for the region (required).
- `-o, --output_dir`: Output directory for saving results (required).
- `--plot_file`: Optional, specify the name of the output plot file (default: `read_counts_plot.png`).

### Example

```bash
exscope -b sample.bam -g Homo_sapiens.GRCh38.104.gff3 -r chr1:150000-160000 -tid ENST00000367770 -o results/
```

This command will:
1. Extract read counts from the BAM file for the specified region.
2. Normalize the read counts.
3. Generate a stacked area plot visualizing read counts across exons.
4. Save the output plot and read counts to the specified output directory.

## Output

ExScope generates the following outputs:
- **Read Counts Plot**: A PNG file visualizing read counts across the specified region.
- **Read Counts Text File**: A text file listing read counts per position.

## GFF3/GTF File Download

If you don’t have a GFF3 or GTF file, you can download one for your species and genome build from:

- **Ensembl FTP**: [Ensembl FTP](https://ftp.ensembl.org/pub/current_gff3/)
- **UCSC Genome Browser**: [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/)

Ensure the file matches your reference genome.

## Troubleshooting

If you encounter any issues, consider the following:
- Ensure that all input files (BAM, GFF3/GTF) are correctly formatted and correspond to the same reference genome.
- Verify that the Ensembl transcript ID matches the specified genomic region.
- Check the log output for warnings or errors.

For further assistance, raise an issue on the [GitHub repository](https://github.com/Taimoor-Khan-bt/exscope).

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request. For major changes, open an issue to discuss your proposed modifications.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
```
