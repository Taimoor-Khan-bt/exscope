import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import logging
import argparse

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract read counts from a BAM file in a specific genomic region.")
    parser.add_argument('-b', '--bam', required=True, help="Input BAM file")
    parser.add_argument('-g', '--gff3', required=True, help="GFF3 file for exon annotations")
    parser.add_argument('-r', '--region', required=True, help="Genomic region in the format 'chr:start-end'")
    parser.add_argument('-tid', '--transcript_id', required=True, help="Transcript id")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory for results")
    parser.add_argument('--plot_file', default="read_counts_plot.png", help="Name of the output plot file (default: read_counts_plot.png)")
    return parser.parse_args()

def read_gff3(gff3_file, region, transcript_id):
    """Read GFF3 file and extract exons from the specified transcript within a given region."""
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))

    try:
        # Read GFF3 file
        gff_df = pd.read_csv(gff3_file, sep="\t", comment='#', header=None,
                             names=["seqname", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])

        # Filter for exons in the specified region
        exons = gff_df[(gff_df["seqname"] == chrom) &
                       (gff_df["feature"] == "exon") &
                       (gff_df["start"] >= start) &
                       (gff_df["end"] <= end)]

        if exons.empty:
            logging.warning(f"No exons found in the specified region: {region}")
            return pd.DataFrame()

        # Parse attributes to get the transcript_id
        def get_transcript_id(attributes):
            """Extract transcript_id from the attributes column."""
            for attr in attributes.split(";"):
                if attr.startswith("transcript_id="):
                    return attr.split("=")[1].strip('"')
            return None

        # Add a column for transcript_id
        exons["transcript_id"] = exons["attributes"].apply(get_transcript_id)

        # Filter for the specified transcript_id
        exons = exons[exons["transcript_id"] == transcript_id]
        
        if exons.empty:
            logging.warning(f"No exons found for transcript ID: {transcript_id}")
        else:
            logging.info(f"Selected transcript: {transcript_id} ({len(exons)} exons)")

        return exons
    
    except Exception as e:
        logging.error(f"Error reading GFF3 file: {e}")
        return pd.DataFrame()

def extract_read_counts(bam_file, region):
    """Extract read counts per base in the specified region."""
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    
    # Open the BAM file
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        logging.error(f"Error opening BAM file: {e}")
        return [], chrom, start, end
    
    # Collect read counts per base in the region
    read_counts = [0] * (end - start + 1)
    for pileupcolumn in bam.pileup(chrom, start, end):
        if start <= pileupcolumn.pos <= end:
            read_counts[pileupcolumn.pos - start] = pileupcolumn.n
    
    bam.close()
    return read_counts, chrom, start, end

def write_read_counts_to_file(read_counts, exons, output_file, chrom, start):
    """Write read counts to a file, grouped by exons if available."""
    with open(output_file, 'w') as f:
        f.write("Position\tReadCount\tExon\n")
        for i, count in enumerate(read_counts):
            position = start + i
            exon = exons[(exons["start"] <= position) & (exons["end"] >= position)]
            exon_label = exon["attributes"].iloc[0] if not exon.empty else "N/A"
            f.write(f"{chrom}:{position}\t{count}\t{exon_label}\n")

    logging.info(f"Read counts written to {output_file}")

def plot_stacked_area_read_counts(read_counts, positions, exon_regions=None, output_plot=None):
    """
    Create a stacked area plot where each exon has a different colored peak.
    """
    fig, ax = plt.subplots(figsize=(14, 7))

    # Ensure read_counts is a 2D array
    read_counts = np.array(read_counts)
    if len(read_counts.shape) == 1:  # If read_counts is 1D, reshape to 2D for consistency
        read_counts = read_counts.reshape(1, -1)

    # Convert positions to a NumPy array
    positions = np.array(positions)

    # Define exon regions with different colors
    if exon_regions:
        exon_colors = plt.colormaps['tab10']  # Use the updated method for colormaps
        # Get enough distinct colors from the colormap
        colors = [exon_colors(i / len(exon_regions)) for i in range(len(exon_regions))]

        # Initialize a zero array for stacking
        cumulative_read_counts = np.zeros_like(read_counts[0])

        # Plot each exon as a stacked area
        for i, (exon_start, exon_end) in enumerate(exon_regions):
            # Extract the indices corresponding to this exon
            exon_mask = (positions >= exon_start) & (positions <= exon_end)
            exon_read_counts = read_counts[:, exon_mask].sum(axis=0)

            # Plot the stacked area for this exon
            ax.fill_between(positions[exon_mask], cumulative_read_counts[exon_mask], 
                            cumulative_read_counts[exon_mask] + exon_read_counts, 
                            color=colors[i], label=f"Exon {i+1}")

            # Update the cumulative read counts
            cumulative_read_counts[exon_mask] += exon_read_counts

    # Set axis labels and title
    ax.set_xlabel('Genomic Position', fontsize=14)
    ax.set_ylabel('Cumulative Read Counts', fontsize=14)
    plt.title("Stacked Area Plot of Read Counts Across Exons", fontsize=16)

    # Customize x-axis labels for readability
    ax.set_xticks(np.linspace(min(positions), max(positions), num=6))
    ax.set_xticklabels([f"{int(x):,}" for x in np.linspace(min(positions), max(positions), num=6)], rotation=0, ha="right", fontsize=10)

    # Add gridlines for better readability
    ax.grid(True, linestyle='--', alpha=0.6)

    # Move legend outside the plot
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=12)

    # Adjust layout to prevent clipping of labels
    fig.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust right side to make space for the legend

    # Save or show plot
    if output_plot:
        plt.savefig(output_plot, bbox_inches='tight')
        print(f"Plot saved to {output_plot}")
    else:
        plt.show()


def create_output_directory(output_dir):
    """Create the output directory if it doesn't exist."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory: {output_dir}")

def main():
    # Parse arguments
    args = parse_arguments()
    
    # Create output directory
    create_output_directory(args.output_dir)
    
    # Extract exons from GFF3 file
    exons = read_gff3(args.gff3, args.region, args.transcript_id)
    
    if exons.empty:
        logging.error("No exons found. Exiting.")
        return
    
    # Extract read counts
    read_counts, chrom, start, end = extract_read_counts(args.bam, args.region)
    
    if not read_counts:
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
