import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import logging
import argparse
from concurrent.futures import ThreadPoolExecutor

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
        # Read GFF3 file efficiently with chunking
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

        # Use .loc[] to avoid SettingWithCopyWarning
        exons.loc[:, "transcript_id"] = exons["attributes"].apply(get_transcript_id)

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
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Collect read counts per base in the region
            read_counts = np.zeros(end - start + 1, dtype=int)
            for pileupcolumn in bam.pileup(chrom, start, end):
                if start <= pileupcolumn.pos <= end:
                    read_counts[pileupcolumn.pos - start] = pileupcolumn.n
    except Exception as e:
        logging.error(f"Error opening BAM file: {e}")
        return None, chrom, start, end
    
    # Log1p transformation and normalization
    read_counts = np.log1p(read_counts)
    min_read, max_read = np.min(read_counts), np.max(read_counts)
    read_counts = 4 * (read_counts - min_read) / (max_read - min_read) - 2
    
    return read_counts, chrom, start, end

def write_read_counts_to_file(read_counts, exons, output_file, chrom, start):
    """Write read counts to a file, grouped by exons if available."""
    with open(output_file, 'w') as f:
        f.write("Position\tReadCount\tExon\n")
        for i, count in enumerate(read_counts):
            position = start + i
            exon = exons[(exons["start"] <= position) & (exons["end"] >= position)]
            exon_label = exon["attributes"].iloc[0] if not exon.empty else "N/A"
            f.write(f"{chrom}:{position}\t{count:.3f}\t{exon_label}\n")

    logging.info(f"Read counts written to {output_file}")

def merge_adjacent_exons(exon_list):
    """
    Merge adjacent exons with the same status (Dup/Del) into a range (e.g., Exon 1-10: Dup).
    """
    merged = []
    start, end, status = exon_list[0]

    for i in range(1, len(exon_list)):
        curr_start, curr_end, curr_status = exon_list[i]
        if curr_status == status and curr_start == end + 1:
            # Merge adjacent exons
            end = curr_end
        else:
            # Save the previous range and start a new one
            if start == end:
                merged.append(f"Exon {start}: {status}")
            else:
                merged.append(f"Exon {start}-{end}: {status}")
            start, end, status = curr_start, curr_end, curr_status

    # Add the last exon or range
    if start == end:
        merged.append(f"Exon {start}: {status}")
    else:
        merged.append(f"Exon {start}-{end}: {status}")

    return merged

def plot_stacked_area_read_counts(read_counts, positions, exon_regions=None, output_plot=None):
    """
    Create a high-quality stacked area plot where each exon has a different colored peak.
    Highlights regions where log1p read counts indicate duplications or deletions.
    Annotates those regions with exon numbers, and adjusts positions to avoid overlaps.
    Widen the main plot and reduce legend text size for better clarity.
    """
    # Increase the figure size to make the plot wider and more readable
    fig, ax = plt.subplots(figsize=(18, 7))  # Increased width for better visibility

    # Convert positions to a NumPy array
    positions = np.array(positions)

    # Plot each exon and highlight duplications/deletions
    cumulative_read_counts = np.zeros_like(read_counts)
    exon_annotations = []
    exon_labels = []
    annotation_offsets = {}

    # Track deleted/duplicated exons for merging adjacent ones
    exon_statuses = []

    # Iterate over exon regions and plot
    for exon_num, (exon_start, exon_end) in enumerate(exon_regions, 1):
        exon_mask = (positions >= exon_start) & (positions <= exon_end)
        exon_read_counts = read_counts[exon_mask]
        is_duplication = False
        is_deletion = False
        
        # Highlight duplications
        duplication_mask = exon_read_counts >= 1.8
        if np.any(duplication_mask):
            ax.fill_between(positions[exon_mask][duplication_mask], cumulative_read_counts[exon_mask][duplication_mask], 
                            cumulative_read_counts[exon_mask][duplication_mask] + exon_read_counts[duplication_mask], 
                            color='red', alpha=0.5)
            is_duplication = True
            logging.info(f"Duplicated region detected: Exon {exon_num} ({exon_start}-{exon_end})")

        # Highlight deletions
        deletion_mask = exon_read_counts <= -1.5
        if np.any(deletion_mask):
            ax.fill_between(positions[exon_mask][deletion_mask], cumulative_read_counts[exon_mask][deletion_mask], 
                            cumulative_read_counts[exon_mask][deletion_mask] + exon_read_counts[deletion_mask], 
                            color='blue', alpha=0.5)
            is_deletion = True
            logging.info(f"Deleted region detected: Exon {exon_num} ({exon_start}-{exon_end})")

        # Add exon status to the list for merging adjacent ones
        if is_duplication:
            exon_statuses.append((exon_num, exon_num, 'Dup'))
        elif is_deletion:
            exon_statuses.append((exon_num, exon_num, 'Del'))

        # Stack regular exon read counts
        ax.fill_between(positions[exon_mask], cumulative_read_counts[exon_mask], 
                        cumulative_read_counts[exon_mask] + exon_read_counts, color='green', alpha=0.3)
        cumulative_read_counts[exon_mask] += exon_read_counts

        # Store exon label for the legend
        exon_labels.append(f"Exon {exon_num}")  # Shortened label for clarity

    # Merge adjacent duplicated/deleted exons and prepare for legend
    if exon_statuses:
        exon_annotations.extend(merge_adjacent_exons(exon_statuses))

    # Set axis labels and title with appropriate fontsize
    ax.set_xlabel('Genomic Position', fontsize=16)
    ax.set_ylabel('Log1p Transformed Read Counts', fontsize=16)
    plt.title(f"Read Count Plot", fontsize=18)

    # Customize x-axis labels for readability
    ax.set_xticks(np.linspace(min(positions), max(positions), num=6))
    ax.set_xticklabels([f"{int(x):,}" for x in np.linspace(min(positions), max(positions), num=6)], rotation=0, ha="right", fontsize=12)

    # Add gridlines for better readability
    ax.grid(True, linestyle='--', alpha=0.6)

    # Build custom legend with smaller font size for clarity
    handles, labels = ax.get_legend_handles_labels()

    # Add exon annotations for duplication and deletion regions
    if exon_annotations:
        for annotation in exon_annotations:
            handles.append(plt.Line2D([0], [0], color='black', lw=4, alpha=0.5))
            labels.append(annotation)

    # Add normal exon regions to the legend
    for label in exon_labels:
        handles.append(plt.Line2D([0], [0], color='green', lw=4, alpha=0.3))
        labels.append(label)

    # Adjust the number of legend columns dynamically based on the number of exons
    ncol = max(1, min(len(exon_labels) // 10 + 1, 5))  # Max 5 columns, adjust based on the number of exons
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1), fontsize=8, ncol=ncol)  # Reduced fontsize to 8 for better clarity

    # Adjust layout to prevent clipping of labels
    fig.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust right side to make space for the legend

    # Save or show plot with a high-quality DPI setting
    if output_plot:
        plt.savefig(output_plot, bbox_inches='tight', dpi=300)  # Set DPI to 300 for high-quality output
        logging.info(f"Plot saved to {output_plot}")
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
    
    if read_counts is None:
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
