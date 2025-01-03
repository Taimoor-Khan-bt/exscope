import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import logging
import argparse
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract read counts from a BAM file in a specific genomic region.")
    parser.add_argument('-b', '--bam', required=True, help="Input BAM file")
    parser.add_argument('-g', '--gff3', required=True, help="GFF3 file for exon annotations")
    parser.add_argument('-r', '--region', required=True, help="Genomic region in the format 'chr:start-end'")
    parser.add_argument('-tid', '--transcript_id', required=True, help="Transcript id")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory for results")
    parser.add_argument('--plot_file', default="read_counts_plot.png", help="Name of the output plot file (default: read_counts_plot.png)")
    return parser.parse_args()

def create_output_directory(output_dir):
    """Create the output directory if it doesn't exist."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory: {output_dir}")

def read_gff3(gff3_file, region, transcript_id):
    """Read and filter the GFF3 file for the specified transcript ID and region."""
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))

    logging.info(f"Reading GFF3 file: {gff3_file}")
    try:
        gff_df = pd.read_csv(gff3_file, sep="\t", comment='#', header=None,
                             names=["seqname", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])
        logging.info(f"Loaded GFF3 file with {gff_df.shape[0]} rows")

        # Standardize chromosome names
        chrom = chrom.replace("chr", "")
        gff_df["seqname"] = gff_df["seqname"].str.replace("chr", "")

        # Filter for overlapping exons
        logging.info("Filtering for overlapping exons in the specified region...")
        exons = gff_df[(gff_df["seqname"] == chrom) &
                       (gff_df["end"] >= start) &
                       (gff_df["start"] <= end) &
                       (gff_df["feature"] == "exon")]

        # Parse transcript ID
        def get_transcript_id(attributes):
            for attr in attributes.split(";"):
                if attr.startswith("transcript_id="):
                    return attr.split("=")[1].strip('"')
            return None

        exons["transcript_id"] = exons["attributes"].apply(get_transcript_id)
        exons_filtered = exons[exons["transcript_id"] == transcript_id]

        if exons_filtered.empty:
            logging.warning(f"No exons found for transcript ID: {transcript_id}")
        else:
            logging.info(f"Selected transcript: {transcript_id} ({len(exons_filtered)} exons)")

        return exons_filtered

    except Exception as e:
        logging.error(f"Error reading GFF3 file: {e}")
        return pd.DataFrame()

def extract_read_counts(bam_file, region):
    """Extract read counts per base in the specified region."""
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    
    logging.info(f"Extracting read counts from BAM file: {bam_file} for region {region}")
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Collect read counts per base in the region
            read_counts = np.zeros(end - start + 1, dtype=int)
            for pileupcolumn in bam.pileup(chrom, start, end):
                if start <= pileupcolumn.pos <= end:
                    read_counts[pileupcolumn.pos - start] = pileupcolumn.n
    except Exception as e:
        logging.error(f"Error extracting read counts: {e}")
        return np.zeros(end - start + 1, dtype=int), chrom, start, end

    return read_counts, chrom, start, end

def write_read_counts_to_file(read_counts, exons_filtered, output_file, chrom, start):
    """Write read counts to a file, grouped by exons if available."""
    logging.info(f"Writing read counts to file: {output_file}")
    with open(output_file, 'w') as f:
        f.write("Position\tReadCount\tExon\n")
        for i, count in tqdm(enumerate(read_counts), desc="Writing read counts to file", total=len(read_counts)):
            position = start + i
            exon = exons_filtered[(exons_filtered["start"] <= position) & (exons_filtered["end"] >= position)]
            try:
                exon_label = exon["attributes"].iloc[0] if not exon.empty else "N/A"
            except Exception as e:
                logging.error(f"Error extracting exon label: {e}")
                exon_label = "N/A"
            f.write(f"{position}\t{count}\t{exon_label}\n")

    logging.info(f"Read counts successfully written to {output_file}")

def normalize_read_counts(read_counts):
    """
    Log1p transform and normalize read counts to the range -2 to +2.
    """
    # Log1p transformation
    logging.info("Log1p transformation applied.")
    read_counts = np.log1p(read_counts)

    # Clip outliers and normalize
    percentile_5 = np.percentile(read_counts, 5)
    percentile_95 = np.percentile(read_counts, 95)
    read_counts = np.clip(read_counts, percentile_5, percentile_95)
    min_read, max_read = np.nanmin(read_counts), np.nanmax(read_counts)
    read_counts = 4 * (read_counts - min_read) / (max_read - min_read) - 2
    logging.info("Data normalized to the range -2 to +2.")

    # Handle NaN values
    read_counts = np.nan_to_num(read_counts, nan=0)  # Replace NaNs with -2 (or any other appropriate value)
    
    return read_counts

def plot_stacked_area(read_counts, positions, exon_regions=None, ax=None):
    """
    Create a stacked area plot for genomic data.
    Highlights duplication (red), deletion (blue), and exons with unique colors.
    """
    if exon_regions is None:
        exon_regions = []

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 8))

    # Ensure positions is a NumPy array
    positions = np.array(positions)

    # Initialize cumulative counts
    cumulative_read_counts = np.zeros_like(read_counts)

    # Define a colormap for exons
    colors = plt.cm.tab20.colors  # A colormap with 20 distinct colors

    # Plot exons with unique colors
    for exon_num, (exon_start, exon_end) in enumerate(exon_regions, 1):
        exon_mask = (positions >= exon_start) & (positions <= exon_end)
        exon_read_counts = read_counts[exon_mask]

        ax.fill_between(
            positions[exon_mask],
            cumulative_read_counts[exon_mask],
            cumulative_read_counts[exon_mask] + exon_read_counts,
            color=colors[exon_num % len(colors)],
            alpha=0.8,  # Increased opacity for better visibility
            linewidth=2,  # Thicker border around the filled areas
            label=f"Exon {exon_num}"
        )
        cumulative_read_counts[exon_mask] += exon_read_counts

    # Aesthetics
    ax.set_xlabel('Genomic Position', fontsize=14)
    ax.set_ylabel('Normalized Log1p Transformed Read Counts', fontsize=14)
    ax.set_title('Stacked Area Plot of Read Counts', fontsize=16)
    ax.set_xticks(np.linspace(min(positions), max(positions), num=8))
    ax.set_xticklabels([f"{int(x):,}" for x in np.linspace(min(positions), max(positions), num=8)],
                       rotation=45, ha='right', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)

    # Add legend outside the plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, title="Exon Regions", ncol=1)

    return ax

def plot_stacked_area_read_counts(read_counts, positions, exon_regions=None, output_plot=None):
    """
    Generate a stacked area plot for the genomic data.
    """
    # Normalize read counts
    read_counts = normalize_read_counts(read_counts)

    # Create subplots
    fig, ax = plt.subplots(figsize=(20, 8))

    # Generate stacked area plot
    plot_stacked_area(read_counts, positions, exon_regions=exon_regions, ax=ax)

    # Adjust layout to fit the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for the legend

    # Save or show the plot
    if output_plot:
        plt.savefig(output_plot, dpi=300, bbox_inches='tight')
        logging.info(f"Plot saved to {output_plot}")
    else:
        plt.show()

def main():
    """
    Main function to parse arguments, create output directory, extract exons and read counts,
    write read counts to a file, and generate a plot.
    """
    # Parse arguments
    args = parse_arguments()
    
    # Create output directory
    create_output_directory(args.output_dir)
    
    # Extract exons from GFF3 file
    exons = read_gff3(args.gff3, args.region, args.transcript_id)
    
    if exons is None or exons.empty:
        logging.error("No exons found. Exiting.")
        return
    
    # Extract read counts
    read_counts, chrom, start, end = extract_read_counts(args.bam, args.region)
    
    if read_counts is None or len(read_counts) == 0:
        logging.error("No read counts extracted. Exiting.")
        return
    
    # Write read counts to file
    output_file = os.path.join(args.output_dir, "read_counts.txt")
    write_read_counts_to_file(read_counts, exons, output_file, chrom, start)
    
    # Prepare exon regions for plotting
    exon_regions = list(zip(exons["start"], exons["end"]))
    
    # Generate the plot
    output_plot = os.path.join(args.output_dir, args.plot_file)
    plot_stacked_area_read_counts(read_counts, list(range(start, end + 1)), exon_regions, output_plot)

if __name__ == "__main__":
    main()