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
from exscope.core.chromosome_coverage import ChromosomeCoverageAnalyzer, ChromosomeAnalysisResult
from exscope.visualizers.chromosome_plot import ChromosomePlotter


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


def analyze_chromosome_coverage(args):
    """Analyze coverage across all chromosomes to detect deletions."""
    logging.info(f"Analyzing chromosome-level coverage")
    
    try:
        # Parse BAM files
        if hasattr(args, 'bams') and args.bams:
            bam_files = args.bams
        elif hasattr(args, 'bams_file') and args.bams_file:
            with open(args.bams_file) as f:
                bam_files = [line.strip() for line in f if line.strip()]
        else:
            bam_files = [args.bam]
        
        logging.info(f"Processing {len(bam_files)} BAM file(s)")
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Analyze each BAM file
        results = []
        for bam_file in bam_files:
            bam_path = Path(bam_file)
            sample_id = bam_path.stem.replace('_deduped', '').replace('.bam', '')
            
            logging.info(f"Analyzing: {sample_id}")
            
            # Initialize analyzer
            analyzer = ChromosomeCoverageAnalyzer(
                bam_path=bam_path,
                reference=args.reference,
                threads=args.threads
            )
            
            # Run analysis
            output_prefix = str(output_dir / f"{sample_id}_chromosome")
            result = analyzer.analyze_chromosomes(
                output_prefix=output_prefix,
                sample_id=sample_id
            )
            
            results.append(result)
            
            # Generate text report
            report_path = output_dir / f"{sample_id}_chromosome_report.txt"
            generate_text_report(result, report_path)
            
            # Save summary table
            summary_table = result.get_summary_table()
            table_path = output_dir / f"{sample_id}_chromosome_summary.tsv"
            summary_table.to_csv(table_path, sep='\t', index=False)
            
            logging.info(f"  Karyotype: {result.inferred_karyotype.value}")
            logging.info(f"  Autosomal mean: {result.autosomal_mean:.2f}x")
            if result.deletions:
                logging.warning(f"  Deletions: {', '.join(result.deletions)}")
        
        # Create visualizations with user-specified dimensions and DPI
        plotter = ChromosomePlotter(width=args.width, height=args.height, dpi=args.dpi)
        
        # Get export formats
        export_formats = args.export_format if hasattr(args, 'export_format') else ['html']
        
        if len(results) == 1:
            # Single sample visualization
            result = results[0]
            
            # Generate plots in all requested formats
            generated_files = []
            
            for fmt in export_formats:
                # Always generate the optimized chromosome tracks (log-scaled for visibility)
                track_path = output_dir / f"{result.sample_id}_chromosome_tracks.{fmt}"
                plotter.plot_chromosome_tracks(result, output=str(track_path))
                generated_files.append(('Chromosome tracks (optimized)', track_path))
                
                # Generate additional plots only if --all-plots is specified
                if args.all_plots:
                    # Bar plot - traditional view (less useful for targeted sequencing)
                    bar_path = output_dir / f"{result.sample_id}_chromosome_coverage.{fmt}"
                    plotter.plot_single_sample(result, output=str(bar_path))
                    generated_files.append(('Bar chart (traditional)', bar_path))
                    
                    # Comprehensive report plot (dual panel)
                    report_plot_path = output_dir / f"{result.sample_id}_chromosome_report.{fmt}"
                    plotter.create_report_plot(result, output=str(report_plot_path))
                    generated_files.append(('Report (dual panel)', report_plot_path))
            
            print(f"\n✓ Chromosome analysis complete for {result.sample_id}")
            print(f"  Karyotype: {result.inferred_karyotype.value}")
            print(f"  Autosomal mean coverage: {result.autosomal_mean:.2f}x")
            
            if result.deletions:
                print(f"  ⚠️  DELETIONS DETECTED: {', '.join(result.deletions)}")
            
            print(f"\nGenerated plots ({', '.join(export_formats)} format{'s' if len(export_formats) > 1 else ''}):")
            for desc, path in generated_files:
                print(f"  {desc}: {path}")
            
            print(f"\nAdditional outputs:")
            print(f"  Text report: {output_dir / f'{result.sample_id}_chromosome_report.txt'}")
            print(f"  Summary table: {output_dir / f'{result.sample_id}_chromosome_summary.tsv'}")
            
            if 'png' in export_formats or 'pdf' in export_formats or 'svg' in export_formats:
                print(f"\n📊 Publication-ready figures generated at {args.dpi} DPI")
            
        else:
            # Multi-sample comparison
            comparison_path = output_dir / "chromosome_comparison.html"
            plotter.plot_comparison(results, output=str(comparison_path))
            
            heatmap_path = output_dir / "chromosome_heatmap.html"
            plotter.plot_heatmap(results, output=str(heatmap_path))
            
            # Combined comparison table
            analyzer_temp = ChromosomeCoverageAnalyzer(bam_files[0])
            comparison_table = analyzer_temp.compare_samples(results)
            comparison_table_path = output_dir / "chromosome_comparison.tsv"
            comparison_table.to_csv(comparison_table_path, sep='\t', index=False)
            
            print(f"\n✓ Chromosome analysis complete for {len(results)} samples")
            print(f"\nOutput files:")
            print(f"  Comparison plot: {comparison_path}")
            print(f"  Heatmap: {heatmap_path}")
            print(f"  Comparison table: {comparison_table_path}")
            
            # Summary of findings
            print(f"\nFindings summary:")
            for result in results:
                status = "✓ Normal" if not result.deletions else f"⚠️  Deletions: {', '.join(result.deletions)}"
                print(f"  {result.sample_id}: {result.inferred_karyotype.value} - {status}")
        
    except Exception as e:
        import traceback
        logging.error(f"Chromosome coverage analysis failed: {e}")
        traceback.print_exc()
        sys.exit(1)


def generate_text_report(result, output_path: Path):
    """Generate human-readable text report with ASCII visualization."""
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"CHROMOSOME COVERAGE ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Sample ID: {result.sample_id}\n")
        f.write(f"Inferred Karyotype: {result.inferred_karyotype.value}\n")
        f.write(f"Autosomal Mean Coverage: {result.autosomal_mean:.2f}x\n\n")
        
        if result.deletions:
            f.write("⚠️  DELETIONS DETECTED:\n")
            for chrom in result.deletions:
                f.write(f"   - {chrom}\n")
            f.write("\n")
        
        if result.warnings:
            f.write("Warnings:\n")
            for warning in result.warnings:
                f.write(f"   - {warning}\n")
            f.write("\n")
        
        f.write("-" * 80 + "\n")
        f.write("CHROMOSOME COVERAGE DETAILS\n")
        f.write("-" * 80 + "\n\n")
        
        # Get max coverage for scaling
        max_cov = max(cov.mean_coverage for cov in result.chromosomes.values())
        
        # Print each chromosome
        chromosomes = sorted(result.chromosomes.keys(), 
                           key=ChromosomeAnalysisResult._sort_chromosomes)
        
        f.write(f"{'Chromosome':<12} {'Coverage':>10} {'Ratio':>8} {'Status':<20} {'Visual':>40}\n")
        f.write("-" * 90 + "\n")
        
        for chrom in chromosomes:
            cov = result.chromosomes[chrom]
            
            # Generate ASCII bar
            bar = ChromosomePlotter.generate_ascii_bar(cov.mean_coverage, max_cov, width=40)
            
            # Status symbol
            status_symbol = {
                'normal': '✓',
                'reduced': '⚠',
                'deleted': '✗',
                'elevated': '↑',
                'partial': '◐',
                'not_sequenced': '○'
            }.get(cov.status.value, ' ')
            
            # Add percent covered info for partial chromosomes
            status_str = cov.status.value
            if cov.status.value == 'partial':
                status_str = f"{cov.status.value} ({cov.percent_covered:.1f}%)"
            elif cov.status.value == 'not_sequenced':
                status_str = "not seq."
            
            f.write(f"{chrom:<12} {cov.mean_coverage:>9.2f}x "
                   f"{cov.normalized_ratio:>7.3f} "
                   f"{status_symbol} {status_str:<18} {bar}\n")
        
        f.write("\n")
        f.write("-" * 80 + "\n")
        f.write("INTERPRETATION\n")
        f.write("-" * 80 + "\n\n")
        
        # Provide interpretation
        if result.inferred_karyotype.value.startswith("XY"):
            f.write("Sample appears to be male (XY karyotype)\n")
        elif result.inferred_karyotype.value.startswith("XX"):
            f.write("Sample appears to be female (XX karyotype)\n")
        elif result.inferred_karyotype.value == "X0 (Turner syndrome)":
            f.write("Sample shows X0 pattern consistent with Turner syndrome\n")
        elif "Klinefelter" in result.inferred_karyotype.value:
            f.write("Sample shows XXY pattern consistent with Klinefelter syndrome\n")
        
        if result.deletions:
            f.write(f"\n⚠️  {len(result.deletions)} chromosome(s) show deletion pattern:\n")
            for chrom in result.deletions:
                cov = result.chromosomes[chrom]
                f.write(f"   {chrom}: Coverage is {cov.normalized_ratio:.1%} of autosomal mean\n")
                f.write(f"          This indicates likely deletion or loss of this chromosome\n")
        
        f.write("\n" + "=" * 80 + "\n")


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
  
  # Analyze chromosome-level coverage (detect Y deletion)
  exscope chromosome-coverage sample.bam -d chr_output/
  
  # Compare multiple samples
  exscope chromosome-coverage --bams sample1.bam sample2.bam -d comparison/
  
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
    
    # Chromosome coverage analysis (deletion detection)
    chr_cov_parser = subparsers.add_parser('chromosome-coverage', 
                                           help='Analyze chromosome-level coverage to detect deletions (lightweight)')
    chr_cov_parser.add_argument('bam', nargs='?', help='Input BAM/CRAM file (single sample mode)')
    bam_group = chr_cov_parser.add_mutually_exclusive_group()
    bam_group.add_argument('--bams', nargs='+', help='Multiple BAM files for comparison')
    bam_group.add_argument('--bams-file', help='File containing BAM file paths (one per line)')
    chr_cov_parser.add_argument('-d', '--output-dir', required=True, help='Output directory')
    chr_cov_parser.add_argument('--reference', help='Reference genome FASTA (for CRAM)')
    chr_cov_parser.add_argument('--all-plots', action='store_true', 
                               help='Generate all plot types (tracks, bar chart, and report). Default: tracks only')
    chr_cov_parser.add_argument('--dpi', type=int, default=300, help='Image resolution for publication (default: 300)')
    chr_cov_parser.add_argument('--threads', type=int, default=1, help='Number of threads (default: 1, max: 2)')
    chr_cov_parser.add_argument('--export-format', nargs='+', 
                               choices=['html', 'png', 'pdf', 'svg'],
                               default=['html'],
                               help='Output format(s) for plots (default: html). Use multiple for publication: --export-format html png pdf')
    chr_cov_parser.add_argument('--width', type=int, default=1400, help='Plot width in pixels (default: 1400)')
    chr_cov_parser.add_argument('--height', type=int, default=800, help='Plot height in pixels (default: 800)')
    chr_cov_parser.set_defaults(func=analyze_chromosome_coverage)
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Execute command
    args.func(args)


if __name__ == "__main__":
    main()
