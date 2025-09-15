import argparse
import sys
import os
from pathlib import Path
import gffutils
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
from tqdm import tqdm
import logging

class BamVisualizer:
    def __init__(self, bam_path, gtf_path, transcript_id, output=None, 
                 dpi=300, reference=None):
        logging.info("Initializing BamVisualizer...")
        self.bam_path = Path(bam_path)
        self.gtf_path = Path(gtf_path)
        self.transcript_id = transcript_id
        self.output = Path(output) if output else None
        if self.output and self.output.is_dir():
            self.output = self.output / f"{self.transcript_id}.png"
        self.dpi = dpi
        self.reference = reference
        self.fig = None
        self.db = None
        self.transcript = None
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.exons = []
        self.coverage = None
        self.reads = []

        self.initialize()
        logging.info("BamVisualizer initialized.")

    def initialize(self):
        """Initialize the visualizer by validating files and parsing GTF"""
        logging.info("Validating files and parsing annotation...")
        self.validate_files()
        self.parse_gtf()
        logging.info("Files validated and annotation parsed.")

    def validate_files(self):
        """Check input file existence and format"""
        logging.info("Validating input files...")
        if not self.bam_path.exists():
            raise FileNotFoundError(f"BAM/CRAM file not found: {self.bam_path}")
        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GTF/GFF file not found: {self.gtf_path}")
        logging.info("Input files validated.")

    def parse_gtf(self):
        """Parse GTF/GFF and extract transcript information with database freshness check"""
        logging.info("Parsing GTF/GFF file...")
        db_path = f"{self.gtf_path}.db"
        if os.path.exists(db_path):
            # Check if GTF is newer than the database
            gtf_mtime = os.path.getmtime(self.gtf_path)
            db_mtime = os.path.getmtime(db_path)
            if gtf_mtime > db_mtime:
                logging.info("GTF file is newer than the database. Removing old database.")
                os.remove(db_path)
        
        try:
            if not os.path.exists(db_path):
                logging.info("Creating new gffutils database...")
                self.db = gffutils.create_db(
                    str(self.gtf_path),
                    dbfn=db_path,
                    force=True,
                    keep_order=True,
                    merge_strategy="create_unique"
                )
                logging.info("Database created.")
            else:
                logging.info("Loading existing gffutils database.")
                self.db = gffutils.FeatureDB(db_path)
            
            logging.info(f"Fetching transcript: {self.transcript_id}")
            self.transcript = self.db[self.transcript_id]
            self.chrom = self.transcript.seqid
            self.start = self.transcript.start
            self.end = self.transcript.end
            self.strand = self.transcript.strand
            self.exons = list(self.db.children(self.transcript, featuretype='exon'))
            if not self.exons:
                raise ValueError(f"No exons found for transcript {self.transcript_id}")
            logging.info(f"Transcript {self.transcript_id} found on chromosome {self.chrom}:{self.start}-{self.end}")
        except gffutils.FeatureNotFoundError:
            raise ValueError(f"Transcript ID {self.transcript_id} not found in GTF")
        logging.info("GTF/GFF file parsed successfully.")

    def load_alignments(self):
        """Load reads from BAM/CRAM file using pysam with region fetch"""
        logging.info("Loading alignments...")
        # Determine file type and open with appropriate parameters
        file_type = 'rc' if self.bam_path.suffix.lower() in ('.cram',) else 'rb'
        open_kwargs = {}
        if file_type == 'rc' and self.reference:
            open_kwargs['reference_filename'] = self.reference
        
        try:
            bam = pysam.AlignmentFile(str(self.bam_path), file_type, **open_kwargs)
        except Exception as e:
            raise RuntimeError(f"Error opening alignment file: {str(e)}")
        
        try:
            if not bam.has_index():
                raise RuntimeError("BAM/CRAM file must be indexed for efficient querying")
            
            # Get total reads in region for progress bar
            total = bam.count(self.chrom, self.start, self.end)
            with tqdm(total=total, desc="Loading reads", unit="reads") as pbar:
                for read in bam.fetch(self.chrom, self.start, self.end):
                    self.reads.append(read)
                    pbar.update(1)
        
            # Sort reads by their start position for consistent visualization
            self.reads.sort(key=lambda x: x.reference_start)
        
        except Exception as e:
            raise RuntimeError(f"Error reading alignment file: {str(e)}")
        finally:
            bam.close()
        logging.info(f"Loaded {len(self.reads)} alignments.")

    def calculate_coverage(self):
        """Calculate strand-specific coverage"""
        logging.info("Calculating coverage...")
        self.coverage = {
            '+': np.zeros(self.end - self.start + 1),
            '-': np.zeros(self.end - self.start + 1)
        }
        
        for read in self.reads:
            if read.is_unmapped:
                continue
            
            strand = '-' if read.is_reverse else '+'
            for pos in (pos for pos in read.get_reference_positions() if pos is not None):
                if self.start <= pos <= self.end:
                    self.coverage[strand][pos - self.start] += 1
        logging.info("Coverage calculated.")

    def plot_genome_view(self):
        """Create multi-track visualization with improved layout"""
        logging.info("Plotting genome view...")
        plt.close('all')  # Close all existing figures
        self.fig = plt.figure(figsize=(15, 10), dpi=self.dpi)
        gs = self.fig.add_gridspec(3, 1, height_ratios=[1, 0.5, 0.5])
        
        ax_cov = self.fig.add_subplot(gs[0])
        ax_cov.fill_between(
            range(len(self.coverage['+'])),
            self.coverage['+'],
            label='Forward',
            color='blue',
            alpha=0.4
        )
        ax_cov.fill_between(
            range(len(self.coverage['-'])) ,
            -self.coverage['-'],  # Negate coverage for reverse strand to plot it below the x-axis
            label='Reverse',
            color='red',
            alpha=0.4
        )
        ax_cov.axhline(y=0, color='black', linewidth=0.8)
        ax_cov.set_title(f"Stranded Coverage - {self.transcript_id} ({self.chrom}:{self.start}-{self.end})")
        ax_cov.set_ylabel("Depth")
        ax_cov.legend()
        ax_cov.set_xlim(0, self.end - self.start)
        
        ax_gene = self.fig.add_subplot(gs[1])
        self.plot_gene_structure(ax_gene)
        
        plt.tight_layout()
        logging.info("Genome view plotted.")

    def plot_gene_structure(self, ax):
        """Visualize exons with labels"""
        height = 0.4
        y_center = 0.5
        
        # Sort exons in transcript order
        sorted_exons = sorted(self.exons, key=lambda x: x.start, reverse=(self.strand == '-'))
        
        # Draw exons with labels
        for i, exon in enumerate(sorted_exons, 1):
            # Calculate exon position relative to transcript start
            exon_start = exon.start - self.start
            exon_end = exon.end - self.start
            ax.add_patch(plt.Rectangle(
                (exon_start, y_center - height/2),
                exon_end - exon_start,
                height,
                facecolor='darkgreen',
                edgecolor='black',
                zorder=2
            ))
            # Add exon number label
            midpoint_x = exon_start + (exon_end - exon_start) / 2
            ax.text(
                midpoint_x, y_center,
                f"Exon {i}",
                ha='center', va='top',
                rotation=45,
                color='black',
                fontsize=6,
                clip_on=True
            )
        
        # Draw introns with direction indicator
        exons_sorted = sorted(self.exons, key=lambda x: x.start)
        for i in range(len(exons_sorted) - 1):
            start = exons_sorted[i].end
            end = exons_sorted[i+1].start
            mid = (start + end) // 2 - self.start
            
            ax.plot(
                [start - self.start, end - self.start],
                [y_center, y_center],
                color='black',
                linestyle='--',
                zorder=1
            )
            
            if self.strand == '+':
                ax.annotate('', xy=(mid + 5, y_center), xytext=(mid -5, y_center),
                            arrowprops=dict(arrowstyle="->", color='black'))
            else:
                ax.annotate('', xy=(mid -5, y_center), xytext=(mid +5, y_center),
                            arrowprops=dict(arrowstyle="->", color='black'))
        
        ax.set_title("Gene Structure with Exon Labels")
        ax.set_xlim(0, self.end - self.start)
        ax.set_ylim(0, 1)
        ax.axis('off')

    def save_or_show(self):
        """Save the figure if an output path is provided; otherwise, inform the user."""
        logging.info("Saving or showing plot...")
        if self.output:
            plt.savefig(self.output, bbox_inches='tight', dpi=self.dpi)
            print(f"Saved visualization to {self.output}")
        else:
            print("No output file specified. To save the plot, please provide an output file path using the -o/--output argument.", file=sys.stderr)
        plt.close(self.fig)
        logging.info("Plot saved or shown.")


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    parser = argparse.ArgumentParser(
        description="Interactive Genomic Visualization Tool with Multi-Panel Layout",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.2.0")
    parser.add_argument("bam", help="Input BAM/CRAM file")
    parser.add_argument("gtf", help="Annotation GTF/GFF file")
    parser.add_argument("-t", "--transcript", required=True, help="Transcript ID to visualize")
    parser.add_argument("-o", "--output", help="Output file (PNG/SVG/PDF)")
    parser.add_argument("--dpi", type=int, default=300, help="Output image resolution")
    parser.add_argument("--reference", help="Reference genome FASTA (required for CRAM)")
    
    args = parser.parse_args()
    
    print(f"BAM: {args.bam}")
    print(f"GTF: {args.gtf}")
    print(f"Transcript: {args.transcript}")
    print(f"Output: {args.output}")
    print(f"DPI: {args.dpi}")
    print(f"Reference: {args.reference}")
    
    try:
        logging.info("Starting visualization process.")
        visualizer = BamVisualizer(
            args.bam,
            args.gtf,
            args.transcript,
            output=args.output,
            dpi=args.dpi,
            reference=args.reference
        )
        visualizer.load_alignments()
        if not visualizer.reads:
            raise ValueError("No reads found in the specified region")
        visualizer.calculate_coverage()
        visualizer.plot_genome_view()
        visualizer.save_or_show()
        logging.info("Visualization process finished successfully.")
    except Exception as e:
        logging.error("An error occurred during visualization.", exc_info=True)
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()