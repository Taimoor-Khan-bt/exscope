#!/usr/bin/env python3
"""Example: Chromosome coverage analysis for Y deletion detection.

This script demonstrates how to use ExScope's chromosome coverage analysis
to detect Y chromosome deletion in WES data. It shows both programmatic
usage and CLI examples.
"""

import sys
from pathlib import Path

# Add parent directory to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent))

from exscope.core.chromosome_coverage import ChromosomeCoverageAnalyzer
from exscope.visualizers.chromosome_plot import ChromosomePlotter


def example_single_sample_analysis():
    """Example: Analyze a single sample for chromosome deletions."""
    print("=" * 80)
    print("EXAMPLE 1: Single Sample Chromosome Analysis")
    print("=" * 80)
    print()
    
    # File paths (replace with your actual files)
    bam_file = "example_patient.bam"
    output_dir = "chr_analysis"
    
    print(f"Input BAM: {bam_file}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Create analyzer
    print("Step 1: Initialize ChromosomeCoverageAnalyzer...")
    analyzer = ChromosomeCoverageAnalyzer(
        bam_path=bam_file,
        threads=4
    )
    
    # Run analysis
    print("Step 2: Calculate chromosome coverage...")
    output_prefix = f"{output_dir}/sample_chromosome"
    result = analyzer.analyze_chromosomes(
        output_prefix=output_prefix,
        sample_id="patient_001"
    )
    
    # Display results
    print("\nStep 3: Analysis Results:")
    print(f"  Sample ID: {result.sample_id}")
    print(f"  Inferred Karyotype: {result.inferred_karyotype.value}")
    print(f"  Autosomal Mean Coverage: {result.autosomal_mean:.2f}x")
    
    if result.deletions:
        print(f"\n  ⚠️  DELETIONS DETECTED: {', '.join(result.deletions)}")
    else:
        print("\n  ✓ No deletions detected")
    
    # Generate summary table
    print("\nStep 4: Generate summary table...")
    summary = result.get_summary_table()
    print(summary)
    
    # Create visualizations
    print("\nStep 5: Create visualizations...")
    plotter = ChromosomePlotter(width=1400, height=700)
    
    # Interactive bar plot
    html_plot = f"{output_dir}/coverage_plot.html"
    plotter.plot_single_sample(result, output=html_plot)
    print(f"  Interactive plot saved: {html_plot}")
    
    # Comprehensive report
    report_plot = f"{output_dir}/coverage_report.html"
    plotter.create_report_plot(result, output=report_plot)
    print(f"  Report plot saved: {report_plot}")
    
    print("\n✓ Analysis complete!")
    print()


def example_multi_sample_comparison():
    """Example: Compare multiple samples to identify Y deletion."""
    print("=" * 80)
    print("EXAMPLE 2: Multi-Sample Comparison")
    print("=" * 80)
    print()
    
    # File paths (replace with your actual files)
    normal_bam = "normal_male.bam"
    patient_bam = "patient_y_deletion.bam"
    output_dir = "comparison"
    
    print(f"Normal sample: {normal_bam}")
    print(f"Patient sample: {patient_bam}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Analyze both samples
    print("Step 1: Analyze normal male control...")
    analyzer1 = ChromosomeCoverageAnalyzer(normal_bam, threads=4)
    result1 = analyzer1.analyze_chromosomes(
        output_prefix=f"{output_dir}/normal_chr",
        sample_id="normal_male"
    )
    
    print("Step 2: Analyze patient sample...")
    analyzer2 = ChromosomeCoverageAnalyzer(patient_bam, threads=4)
    result2 = analyzer2.analyze_chromosomes(
        output_prefix=f"{output_dir}/patient_chr",
        sample_id="patient"
    )
    
    # Compare results
    print("\nStep 3: Compare Results:")
    print(f"  Normal: {result1.inferred_karyotype.value} - Deletions: {result1.deletions}")
    print(f"  Patient: {result2.inferred_karyotype.value} - Deletions: {result2.deletions}")
    
    # Create comparison visualizations
    print("\nStep 4: Create comparison plots...")
    plotter = ChromosomePlotter(width=1400, height=700)
    
    # Side-by-side comparison
    comparison_plot = f"{output_dir}/comparison.html"
    plotter.plot_comparison([result1, result2], output=comparison_plot)
    print(f"  Comparison plot saved: {comparison_plot}")
    
    # Heatmap
    heatmap_plot = f"{output_dir}/heatmap.html"
    plotter.plot_heatmap([result1, result2], output=heatmap_plot)
    print(f"  Heatmap saved: {heatmap_plot}")
    
    # Comparison table
    comparison_table = analyzer1.compare_samples([result1, result2])
    print("\nStep 5: Comparison Table:")
    print(comparison_table)
    
    print("\n✓ Comparison complete!")
    print()


def example_cli_commands():
    """Example: Show equivalent CLI commands."""
    print("=" * 80)
    print("EXAMPLE 3: CLI Commands")
    print("=" * 80)
    print()
    
    print("Single sample analysis:")
    print("  exscope chromosome-coverage patient.bam -d chr_analysis/")
    print()
    
    print("Multi-sample comparison:")
    print("  exscope chromosome-coverage --bams normal.bam patient.bam -d comparison/")
    print()
    
    print("Batch analysis from file:")
    print("  # Create bam_files.txt with one BAM path per line")
    print("  exscope chromosome-coverage --bams-file bam_files.txt -d batch_chr/")
    print()
    
    print("With CRAM files:")
    print("  exscope chromosome-coverage sample.cram -d output/ --reference hg38.fa")
    print()


def example_interpreting_results():
    """Example: How to interpret chromosome coverage results."""
    print("=" * 80)
    print("EXAMPLE 4: Interpreting Results")
    print("=" * 80)
    print()
    
    print("Coverage Ratio Interpretation:")
    print("  • Ratio ~1.0 for autosomes (chr1-22): NORMAL - expected 2 copies")
    print("  • Ratio ~0.5 for X in males: NORMAL - expected 1 copy")
    print("  • Ratio ~0.5 for Y in males: NORMAL - expected 1 copy")
    print("  • Ratio ~1.0 for X in females: NORMAL - expected 2 copies")
    print("  • Ratio ~0.0 for Y in females: NORMAL - females have no Y chromosome")
    print()
    
    print("Abnormal Patterns:")
    print("  • Ratio <0.15: DELETED - chromosome likely absent")
    print("  • Ratio 0.15-0.65: REDUCED - partial deletion or low coverage")
    print("  • Ratio >1.45: ELEVATED - possible duplication or trisomy")
    print()
    
    print("Karyotype Interpretations:")
    print("  • XX: Normal female")
    print("  • XY: Normal male")
    print("  • X0: Turner syndrome (monosomy X)")
    print("  • XXY: Klinefelter syndrome")
    print("  • XYY: XYY syndrome")
    print("  • XXX: Triple X syndrome")
    print()
    
    print("Y Chromosome Deletion in WES:")
    print("  If a male patient shows:")
    print("    - Autosomal coverage: ~50x (ratio ~1.0)")
    print("    - X chromosome: ~25x (ratio ~0.5)")
    print("    - Y chromosome: <5x (ratio <0.1)")
    print("  → Likely Y chromosome deletion")
    print()
    
    print("Common WES Coverage Patterns:")
    print("  • Uniform coverage across chr1-22: Good quality WES")
    print("  • MT (mitochondrial) may show very high or variable coverage: NORMAL")
    print("  • Sex chromosomes should follow expected ratios based on sex")
    print()


def main():
    """Run all examples."""
    print("\n" + "=" * 80)
    print("ExScope Chromosome Coverage Analysis Examples")
    print("=" * 80)
    print()
    
    print("This script demonstrates chromosome-level coverage analysis")
    print("for detecting deletions and chromosomal abnormalities.")
    print()
    
    # Show CLI commands first (most practical)
    example_cli_commands()
    
    # Show interpretation guide
    example_interpreting_results()
    
    print("=" * 80)
    print("For programmatic usage examples, see the code in this file:")
    print(f"  {__file__}")
    print()
    print("To run actual analysis, uncomment the following lines and")
    print("provide your own BAM files:")
    print()
    
    # Uncomment these to run actual analysis (requires BAM files)
    # example_single_sample_analysis()
    # example_multi_sample_comparison()
    
    print("=" * 80)
    print()


if __name__ == "__main__":
    main()
