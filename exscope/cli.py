#!/usr/bin/env python3
"""ExScope: Enhanced command-line interface with batch processing and advanced features.

This module provides a user-friendly CLI for genomic visualization with support for:
- Single gene visualization
- Batch processing of multiple genes
- Whole chromosome visualization
- Integration with established tools (mosdepth, samtools, pyGenomeTracks)
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import List

from exscope.core.coverage import CoverageCalculator, SequencingType
from exscope.core.batch_processor import BatchProcessor
from exscope.visualizers.genomics_track import GenomeVisualizer, GenomicRegion
from exscope.core.annotations import AnnotationProcessor


def setup_logging(verbose: bool = False):
    """Configure logging based on verbosity level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def visualize_single_gene(args):
    """Visualize a single gene."""
    logging.info(f"Visualizing gene: {args.gene}")
    
    try:
        # Initialize annotation processor
        annotation_proc = AnnotationProcessor(args.annotation)
        
        # Extract gene region
        gene_region = None
        if hasattr(annotation_proc, 'bedtool') and annotation_proc.bedtool:
            filtered = []
            for entry in annotation_proc.bedtool:
                if entry.name == args.gene or entry.name.startswith(f"{args.gene}_"):
                    filtered.append(entry)
            
            if filtered:
                gene_region = GenomicRegion(
                    chrom=filtered[0].chrom,
                    start=min(int(f.start) for f in filtered),
                    end=max(int(f.end) for f in filtered)
                )
        
        if not gene_region:
            raise ValueError(f"Gene {args.gene} not found in annotation file")
        
        logging.info(f"Gene region: {gene_region}")
        
        # Calculate coverage
        calculator = CoverageCalculator(
            bam_path=args.bam,
            seq_type=SequencingType.PANEL,
            target_regions=args.annotation,
            reference=args.reference,
            threads=args.threads
        )
        
        output_prefix = args.output.replace('.png', '_coverage')
        coverage_stats = calculator.calculate_coverage(output_prefix)
        gene_coverage = calculator.get_gene_coverage(output_prefix, args.gene)
        
        # Create visualization
        visualizer = GenomeVisualizer(dpi=args.dpi)
        
        coverage_file = Path(f"{output_prefix}.regions.bed.gz")
        if coverage_file.exists():
            visualizer.add_coverage_track(
                coverage_file,
                title=f"{args.gene} Coverage ({gene_coverage['coverage'].mean():.1f}x)"
            )
        
        visualizer.add_coverage_track(args.annotation, title="Target Regions", color="#1F78B4")
        
        visualizer.plot(region=gene_region, output=args.output, width=40, height_per_track=2)
        
        print(f"\n✓ Visualization saved to: {args.output}")
        print(f"  Mean coverage: {gene_coverage['coverage'].mean():.1f}x")
        print(f"  Regions: {len(gene_coverage)}")
        
    except Exception as e:
        import traceback
        logging.error(f"Failed to visualize gene: {e}")
        traceback.print_exc()
        sys.exit(1)


def visualize_batch(args):
    """Visualize multiple genes in batch mode."""
    # Parse gene list
    if args.genes_file:
        with open(args.genes_file) as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        genes = [g.strip() for g in args.genes.split(',')]
    
    logging.info(f"Batch processing {len(genes)} genes")
    
    try:
        processor = BatchProcessor(
            bam_path=args.bam,
            annotation_path=args.annotation,
            output_dir=args.output_dir,
            reference=args.reference,
            threads=args.threads,
            dpi=args.dpi
        )
        
        results = processor.process_genes(genes, parallel=not args.no_parallel)
        processor.generate_summary_report(results)
        
        successful = sum(1 for r in results if r['status'] == 'success')
        failed = len(results) - successful
        
        print(f"\n✓ Batch processing complete!")
        print(f"  Successful: {successful}")
        print(f"  Failed: {failed}")
        print(f"  Output directory: {args.output_dir}")
        
    except Exception as e:
        logging.error(f"Batch processing failed: {e}")
        sys.exit(1)


def visualize_chromosome(args):
    """Visualize an entire chromosome."""
    logging.info(f"Visualizing chromosome: {args.chromosome}")
    
    try:
        # For chromosome visualization, we'll use the annotation file to show target regions
        # but calculate coverage across the whole chromosome with larger windows
        calculator = CoverageCalculator(
            bam_path=args.bam,
            seq_type=SequencingType.PANEL,  # Use PANEL to work with target regions
            target_regions=args.annotation,
            reference=args.reference,
            threads=args.threads
        )
        
        output_prefix = args.output.replace('.png', '_coverage')
        coverage_stats = calculator.calculate_coverage(output_prefix)
        
        # Create visualization
        visualizer = GenomeVisualizer(dpi=args.dpi)
        
        coverage_file = Path(f"{output_prefix}.regions.bed.gz")
        if coverage_file.exists():
            visualizer.add_coverage_track(
                coverage_file,
                title=f"{args.chromosome} Coverage ({coverage_stats.get('mean_coverage', 0):.1f}x)"
            )
        
        # Add annotation track to show target regions
        visualizer.add_coverage_track(args.annotation, title="Target Regions", color="#1F78B4")
        
        visualizer.plot_chromosome(args.chromosome, output=args.output, width=60, height_per_track=2)
        
        print(f"\n✓ Chromosome visualization saved to: {args.output}")
        print(f"  Mean coverage: {coverage_stats.get('mean_coverage', 0):.1f}x")
        
    except Exception as e:
        import traceback
        logging.error(f"Failed to visualize chromosome: {e}")
        traceback.print_exc()
        sys.exit(1)


def main():
    """Main entry point for ExScope CLI."""
    parser = argparse.ArgumentParser(
        description="ExScope: Genomic Visualization Tool with Established Algorithms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualize a single gene
  exscope single sample.bam targets.bed -g BRCA1 -o brca1.png
  
  # Batch process multiple genes
  exscope batch sample.bam targets.bed -g "BRCA1,BRCA2,TP53" -d output/
  
  # Visualize entire chromosome
  exscope chromosome sample.bam targets.bed -c chr17 -o chr17.png
  
  # Use genes from file
  exscope batch sample.bam targets.bed --genes-file genes.txt -d output/
        """
    )
    
    parser.add_argument("-v", "--version", action="version", version="ExScope 1.0.0")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    subparsers.required = True
    
    # Single gene visualization
    single_parser = subparsers.add_parser('single', help='Visualize a single gene')
    single_parser.add_argument('bam', help='Input BAM/CRAM file')
    single_parser.add_argument('annotation', help='Annotation file (BED/GTF/GFF)')
    single_parser.add_argument('-g', '--gene', required=True, help='Gene name to visualize')
    single_parser.add_argument('-o', '--output', required=True, help='Output file (PNG/SVG/PDF)')
    single_parser.add_argument('--reference', help='Reference genome FASTA (for CRAM)')
    single_parser.add_argument('--dpi', type=int, default=100, help='Image resolution (default: 100)')
    single_parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    single_parser.set_defaults(func=visualize_single_gene)
    
    # Batch processing
    batch_parser = subparsers.add_parser('batch', help='Process multiple genes in batch')
    batch_parser.add_argument('bam', help='Input BAM/CRAM file')
    batch_parser.add_argument('annotation', help='Annotation file (BED/GTF/GFF)')
    batch_group = batch_parser.add_mutually_exclusive_group(required=True)
    batch_group.add_argument('-g', '--genes', help='Comma-separated list of genes')
    batch_group.add_argument('--genes-file', help='File containing gene names (one per line)')
    batch_parser.add_argument('-d', '--output-dir', required=True, help='Output directory')
    batch_parser.add_argument('--reference', help='Reference genome FASTA (for CRAM)')
    batch_parser.add_argument('--dpi', type=int, default=100, help='Image resolution (default: 100)')
    batch_parser.add_argument('--threads', type=int, default=4, help='Number of parallel processes (default: 4)')
    batch_parser.add_argument('--no-parallel', action='store_true', help='Disable parallel processing')
    batch_parser.set_defaults(func=visualize_batch)
    
    # Chromosome visualization
    chr_parser = subparsers.add_parser('chromosome', help='Visualize entire chromosome')
    chr_parser.add_argument('bam', help='Input BAM/CRAM file')
    chr_parser.add_argument('annotation', help='Annotation file (BED/GTF/GFF)')
    chr_parser.add_argument('-c', '--chromosome', required=True, help='Chromosome name (e.g., chr17)')
    chr_parser.add_argument('-o', '--output', required=True, help='Output file (PNG/SVG/PDF)')
    chr_parser.add_argument('--reference', help='Reference genome FASTA (for CRAM)')
    chr_parser.add_argument('--dpi', type=int, default=100, help='Image resolution (default: 100)')
    chr_parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    chr_parser.set_defaults(func=visualize_chromosome)
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Execute command
    args.func(args)


if __name__ == "__main__":
    main()
