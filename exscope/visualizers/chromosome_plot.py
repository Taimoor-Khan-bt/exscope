"""Plotly-based visualization for chromosome-level coverage analysis.

This module creates interactive, publication-quality visualizations for chromosome
coverage comparisons. Particularly useful for detecting deletions and aneuploidies
in WES data.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import pandas as pd
import numpy as np

from exscope.core.chromosome_coverage import (
    ChromosomeAnalysisResult,
    ChromosomeStatus,
    ChromosomeCoverage
)


class ChromosomePlotter:
    """Create interactive visualizations for chromosome coverage analysis."""
    
    # Color scheme
    COLORS = {
        ChromosomeStatus.NORMAL: '#2ecc71',          # Green - Normal coverage
        ChromosomeStatus.REDUCED: '#f39c12',         # Orange - Reduced coverage
        ChromosomeStatus.DELETED: '#e74c3c',         # Red - Deleted/absent
        ChromosomeStatus.ELEVATED: '#3498db',        # Blue - Elevated coverage
        ChromosomeStatus.PARTIAL: '#9b59b6',         # Purple - Partially sequenced
        ChromosomeStatus.NOT_SEQUENCED: '#95a5a6'    # Gray - Not in BAM
    }
    
    # Expected WES coverage per chromosome (approximate % of chromosome covered by exons)
    # Based on RefSeq exon coverage for hg38
    EXPECTED_WES_COVERAGE = {
        'chr1': 2.1, 'chr2': 2.0, 'chr3': 1.9, 'chr4': 1.5, 'chr5': 1.8,
        'chr6': 1.9, 'chr7': 2.0, 'chr8': 1.7, 'chr9': 2.0, 'chr10': 1.8,
        'chr11': 2.3, 'chr12': 2.1, 'chr13': 1.2, 'chr14': 1.5, 'chr15': 1.6,
        'chr16': 2.2, 'chr17': 2.8, 'chr18': 1.3, 'chr19': 3.1, 'chr20': 2.0,
        'chr21': 1.4, 'chr22': 2.5, 'chrX': 1.6, 'chrY': 0.5, 'chrM': 100.0
    }
    
    def __init__(self, width: int = 1200, height: int = 600, dpi: int = 100):
        """Initialize ChromosomePlotter.
        
        Args:
            width: Plot width in pixels
            height: Plot height in pixels
            dpi: DPI for static image export
        """
        self.width = width
        self.height = height
        self.dpi = dpi
    
    @staticmethod
    def get_ratio_color(ratio: float) -> str:
        """Get color for normalized depth ratio based on biological interpretation.
        
        Color scheme based on clinical significance:
        - Green: Normal diploid (0.9-1.1)
        - Yellow: Haploid/hemizygous (0.4-0.6)
        - Orange: Partial loss/gain (0.7-0.9, 1.1-1.3)
        - Red: Complete deletion (<0.3) or duplication (>1.4)
        - Gray: No coverage (0.0)
        
        Args:
            ratio: Normalized coverage depth ratio
            
        Returns:
            Hex color code
        """
        if ratio == 0.0:
            return '#95a5a6'  # Gray - not sequenced
        elif ratio < 0.3:
            return '#e74c3c'  # Red - complete deletion
        elif 0.4 <= ratio <= 0.6:
            return '#f39c12'  # Yellow/Orange - haploid (expected for chrY in males, hemizygous deletions)
        elif 0.7 <= ratio < 0.9:
            return '#e67e22'  # Orange - partial loss
        elif 0.9 <= ratio <= 1.1:
            return '#2ecc71'  # Green - normal diploid
        elif 1.1 < ratio <= 1.3:
            return '#e67e22'  # Orange - partial gain
        elif ratio > 1.4:
            return '#e74c3c'  # Red - duplication/trisomy
        else:  # 0.3-0.4, 0.6-0.7, 1.3-1.4
            return '#c0392b'  # Dark red - abnormal intermediate values
    
    @staticmethod
    def get_ratio_interpretation(ratio: float, chrom: str, karyotype = "Unknown") -> str:
        """Get biological interpretation hint for a normalized depth ratio.
        
        Args:
            ratio: Normalized coverage depth ratio
            chrom: Chromosome name
            karyotype: Inferred karyotype (e.g., "XX", "XY") - can be string or Karyotype enum
            
        Returns:
            Interpretation string
        """
        # Handle both string and Karyotype enum
        karyotype_str = str(karyotype.value) if hasattr(karyotype, 'value') else str(karyotype)
        
        is_chrY = chrom in ['chrY', 'Y']
        is_chrX = chrom in ['chrX', 'X']
        is_male = 'XY' in karyotype_str
        
        if ratio == 0.0:
            return "Not sequenced"
        elif ratio < 0.3:
            return "Complete deletion (nullisomy)"
        elif 0.4 <= ratio <= 0.6:
            if is_chrY and is_male:
                return "Normal (haploid Y in male)"
            elif is_chrX and is_male:
                return "Normal (haploid X in male)"
            else:
                return "Hemizygous deletion (50% loss)"
        elif 0.7 <= ratio < 0.9:
            return "Partial deletion (mosaic or segmental)"
        elif 0.9 <= ratio <= 1.1:
            return "Normal diploid"
        elif 1.1 < ratio <= 1.3:
            return "Partial duplication (mosaic or segmental)"
        elif 1.4 <= ratio < 1.6:
            return "Trisomy (3 copies)"
        elif ratio >= 1.6:
            return "Tetrasomy or higher (≥4 copies)"
        else:
            return "Abnormal - uncertain significance"
    
    @staticmethod
    def _transform_coverage_for_visibility(percent_covered: float, method: str = 'log') -> float:
        """Transform coverage percentage for better visibility of small values.
        
        Args:
            percent_covered: Actual percentage of chromosome covered (0-100)
            method: Transformation method ('log', 'sqrt', 'linear')
            
        Returns:
            Transformed percentage for visualization (0-100)
        """
        if percent_covered <= 0:
            return 0.0
        
        if method == 'log':
            # Logarithmic scaling: makes small values visible while preserving ordering
            # Formula: log10(1 + percent * scale) / log10(1 + 100 * scale) * 100
            # Scale factor of 10 provides good visibility for 0.01% - 10% range
            scale = 10
            transformed = (np.log10(1 + percent_covered * scale) / 
                          np.log10(1 + 100 * scale)) * 100
            return transformed
        
        elif method == 'sqrt':
            # Square root scaling: less aggressive than log
            return np.sqrt(percent_covered) * 10
        
        else:  # linear
            return percent_covered
    
    def plot_coverage_breadth_qc(self, result: ChromosomeAnalysisResult,
                                 output: Optional[str] = None,
                                 show_plot: bool = False) -> go.Figure:
        """Create coverage breadth QC visualization (secondary/quality control metric).
        
        **WARNING**: This is a QC metric, NOT for clinical interpretation of CNVs/aneuploidies.
        Use plot_normalized_depth_ratio() for clinical CNV detection.
        
        This shows the breadth of coverage (% of chromosome with any reads) normalized
        by expected coverage for each chromosome. Useful for:
        - Verifying target capture efficiency in WES
        - Detecting systematic biases in library prep
        - QC for targeted sequencing panels
        
        Normalized breadth = (actual % covered) / (expected % for this chromosome)
        - Value of 1.0 = achieved expected coverage breadth
        - Value < 0.8 = poor capture efficiency
        - Value > 1.2 = better than expected (unusual for WES)
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            
        Returns:
            Plotly Figure object
        """
        # Sort chromosomes
        sorted_chroms = sorted(result.chromosomes.keys(), 
                             key=ChromosomeAnalysisResult._sort_chromosomes)
        
        # Create figure
        fig = go.Figure()
        
        # Prepare data arrays for efficient plotting
        all_shapes = []  # Collect all shapes to add at once
        
        # Add each chromosome as a track
        for idx, chrom in enumerate(reversed(sorted_chroms)):  # Reverse for top-to-bottom display
            cov = result.chromosomes[chrom]
            y_pos = idx
            
            # Get actual percentage covered
            actual_percent = cov.percent_covered
            
            # Get expected coverage for this chromosome (for WES)
            expected_percent = self.EXPECTED_WES_COVERAGE.get(chrom, 2.0)
            
            # Normalize: actual / expected (no arbitrary transformation)
            if expected_percent > 0:
                normalized_breadth = (actual_percent / expected_percent) * 100
            else:
                normalized_breadth = 0
            
            # Cap at 100% for visualization
            visual_percent = min(normalized_breadth, 100)
            
            # Color based on normalized breadth (QC metric)
            if normalized_breadth >= 80:
                color = '#2ecc71'  # Green - good capture
            elif normalized_breadth >= 50:
                color = '#f39c12'  # Orange - moderate capture
            else:
                color = '#e74c3c'  # Red - poor capture
            
            # Add gray background (100% = expected coverage achieved)
            all_shapes.append(dict(
                type="rect",
                x0=0, x1=100,
                y0=y_pos - 0.4, y1=y_pos + 0.4,
                fillcolor='#ecf0f1',
                line=dict(color='#bdc3c7', width=1),
                layer='below'
            ))
            
            # Add colored fill for normalized breadth
            if visual_percent > 0:
                all_shapes.append(dict(
                    type="rect",
                    x0=0, x1=visual_percent,
                    y0=y_pos - 0.4, y1=y_pos + 0.4,
                    fillcolor=color,
                    line=dict(color=color, width=0),
                    opacity=0.8
                ))
            
            # Add reference line at 100% (expected)
            all_shapes.append(dict(
                type="line",
                x0=100, x1=100,
                y0=y_pos - 0.5, y1=y_pos + 0.5,
                line=dict(color='rgba(0,0,0,0.5)', width=2, dash='solid')
            ))
            
            # Add invisible scatter trace for hover information (one point per chromosome)
            hover_text = (
                f"<b>{chrom}</b><br>"
                f"<b>QC Metric - Breadth of Coverage</b><br>"
                f"Actual: {actual_percent:.3f}% of chromosome<br>"
                f"Expected: {expected_percent:.1f}% (WES)<br>"
                f"Normalized: {normalized_breadth:.1f}% of expected<br>"
                f"Mean depth: {cov.mean_coverage:.2f}x<br>"
                f"Covered bases: {cov.covered_bases:,} / {cov.total_bases:,} bp"
            )
            
            fig.add_trace(go.Scatter(
                x=[visual_percent / 2],  # Center of filled region
                y=[y_pos],
                mode='markers',
                marker=dict(size=0.1, opacity=0),  # Invisible marker
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False
            ))
        
        # Add all shapes at once (much faster than individual additions)
        fig.update_layout(shapes=all_shapes)
        
        # Add chromosome labels as y-axis tick labels
        chrom_labels = [result.chromosomes[c].chrom for c in reversed(sorted_chroms)]
        
        # Update layout with proper margins and axis settings
        fig.update_layout(
            title=dict(
                text=f"<b>QC Metric:</b> Coverage Breadth Analysis - {result.sample_id}<br>"
                     f"<sub>Normalized breadth = (actual % covered) / (expected % for chromosome)<br>"
                     f"<b>⚠️ WARNING:</b> This is a QC metric, NOT for clinical CNV/aneuploidy detection<br>"
                     f"Use Normalized Depth Ratio plot for clinical interpretation</sub>",
                x=0.5,
                xanchor='center'
            ),
            xaxis=dict(
                title="Normalized Breadth (% of expected coverage achieved)",
                showgrid=True,
                gridcolor='#f0f0f0',
                zeroline=False,
                range=[0, 110]  # 0-100% with small padding
            ),
            yaxis=dict(
                showticklabels=True,
                tickmode='array',
                tickvals=list(range(len(sorted_chroms))),
                ticktext=chrom_labels,
                tickfont=dict(size=10),
                showgrid=False,
                zeroline=False,
                range=[-0.5, len(sorted_chroms) - 0.5]
            ),
            height=max(600, len(sorted_chroms) * 25),
            width=self.width,
            plot_bgcolor='white',
            hovermode='closest',
            margin=dict(l=80, r=40, t=160, b=100),
            annotations=[
                dict(
                    text="<span style='color:#2ecc71'>■</span> Good (≥80%)  "
                         "<span style='color:#f39c12'>■</span> Moderate (50-80%)  "
                         "<span style='color:#e74c3c'>■</span> Poor (<50%)",
                    xref="paper", yref="paper",
                    x=0.5, y=-0.08,
                    showarrow=False,
                    font=dict(size=10),
                    xanchor='center'
                )
            ]
        )
        
        # Save if output specified
        if output:
            self._save_figure(fig, output)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def plot_chromosome_tracks(self, result: ChromosomeAnalysisResult,
                              output: Optional[str] = None,
                              show_plot: bool = False,
                              transformation: str = 'log') -> go.Figure:
        """DEPRECATED: Use plot_normalized_depth_ratio() for clinical analysis.
        
        This method is deprecated and kept only for backward compatibility.
        It uses misleading log transformations that create false visual impressions.
        
        For clinical CNV/aneuploidy detection, use:
        - plot_normalized_depth_ratio() - PRIMARY clinical visualization
        - plot_coverage_breadth_qc() - QC metric for capture efficiency
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            transformation: Transform method (deprecated parameter, ignored)
            
        Returns:
            Plotly Figure object
        """
        import warnings
        warnings.warn(
            "plot_chromosome_tracks() is deprecated due to misleading log transformations. "
            "Use plot_normalized_depth_ratio() for clinical analysis or "
            "plot_coverage_breadth_qc() for QC metrics.",
            DeprecationWarning,
            stacklevel=2
        )
        # Redirect to the QC breadth plot
        return self.plot_coverage_breadth_qc(result, output, show_plot)
    
    def plot_normalized_depth_ratio(self, result: ChromosomeAnalysisResult,
                                    output: Optional[str] = None,
                                    show_plot: bool = False) -> go.Figure:
        """Create normalized depth ratio plot - the scientifically validated primary visualization.
        
        This visualization shows chromosome coverage normalized by autosomal mean depth:
        - Y-axis: Normalized ratio (chromosome depth / autosomal mean)
        - X-axis: Chromosomes (chr1-22, X, Y, M)
        - Reference lines: 1.0 (diploid), 0.5 (haploid), 1.5 (trisomy)
        - Color coding based on biological interpretation (not arbitrary status)
        
        This is the clinical standard for CNV detection and aneuploidy detection.
        Ratio interpretation:
        - 1.0 ± 0.1: Normal diploid (2 copies)
        - 0.5 ± 0.1: Haploid (1 copy - normal for chrX/chrY in males, or hemizygous deletion)
        - 0.0-0.3: Complete deletion (nullisomy)
        - 1.5 ± 0.1: Trisomy (3 copies)
        - >1.6: Tetrasomy or higher
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            
        Returns:
            Plotly Figure object
        """
        # Prepare data
        chromosomes = []
        ratios = []
        colors = []
        hover_texts = []
        
        sorted_chroms = sorted(result.chromosomes.keys(),
                              key=ChromosomeAnalysisResult._sort_chromosomes)
        
        for chrom in sorted_chroms:
            cov = result.chromosomes[chrom]
            chromosomes.append(chrom)
            ratios.append(cov.normalized_ratio)
            
            # Get color based on biological interpretation of ratio
            colors.append(self.get_ratio_color(cov.normalized_ratio))
            
            # Get interpretation
            interpretation = self.get_ratio_interpretation(
                cov.normalized_ratio, 
                chrom, 
                result.inferred_karyotype.value
            )
            
            # Create detailed hover text
            hover_text = (
                f"<b>{chrom}</b><br>"
                f"Normalized Ratio: {cov.normalized_ratio:.3f}<br>"
                f"Mean Depth: {cov.mean_coverage:.2f}x<br>"
                f"Autosomal Mean: {result.autosomal_mean:.2f}x<br>"
                f"<b>Interpretation: {interpretation}</b><br>"
                f"Covered Bases: {cov.covered_bases:,} / {cov.total_bases:,}<br>"
                f"Status: {cov.status.value}"
            )
            hover_texts.append(hover_text)
        
        # Create figure
        fig = go.Figure()
        
        # Add bar chart
        fig.add_trace(go.Bar(
            x=chromosomes,
            y=ratios,
            marker_color=colors,
            hovertemplate='%{hovertext}<extra></extra>',
            hovertext=hover_texts,
            showlegend=False
        ))
        
        # Add reference lines with biological meaning
        fig.add_hline(
            y=1.0,
            line_dash="solid",
            line_color="green",
            line_width=2,
            annotation_text="Diploid (2 copies)",
            annotation_position="right",
            annotation=dict(font_size=10, font_color="green")
        )
        
        fig.add_hline(
            y=0.5,
            line_dash="dash",
            line_color="orange",
            line_width=1.5,
            annotation_text="Haploid (1 copy)",
            annotation_position="right",
            annotation=dict(font_size=10, font_color="orange")
        )
        
        fig.add_hline(
            y=1.5,
            line_dash="dot",
            line_color="red",
            line_width=1.5,
            annotation_text="Trisomy (3 copies)",
            annotation_position="right",
            annotation=dict(font_size=10, font_color="red")
        )
        
        # Update layout
        title_text = f"Chromosome Coverage Analysis: {result.sample_id}"
        subtitle_parts = [f"Karyotype: {result.inferred_karyotype.value}"]
        subtitle_parts.append(f"Autosomal Mean Depth: {result.autosomal_mean:.2f}x")
        
        if result.deletions:
            subtitle_parts.append(f"⚠️ Deletions: {', '.join(result.deletions)}")
        
        subtitle = " | ".join(subtitle_parts)
        
        fig.update_layout(
            title=dict(
                text=f"{title_text}<br><sub>{subtitle}</sub>",
                x=0.5,
                xanchor='center'
            ),
            xaxis_title="Chromosome",
            yaxis_title="Normalized Coverage Ratio (Chromosome Depth / Autosomal Mean)",
            width=self.width,
            height=self.height,
            template="plotly_white",
            hovermode='x',
            font=dict(size=12),
            yaxis=dict(
                range=[0, max(2.0, max(ratios) * 1.1)],  # Dynamic range, minimum 0-2
                gridcolor='#f0f0f0'
            ),
            annotations=[
                dict(
                    text="<b>Ratio Interpretation:</b> "
                         "<span style='color:#2ecc71'>■</span> Normal (0.9-1.1)  "
                         "<span style='color:#f39c12'>■</span> Haploid (0.4-0.6)  "
                         "<span style='color:#e67e22'>■</span> Partial Loss/Gain  "
                         "<span style='color:#e74c3c'>■</span> Deletion (<0.3) or Duplication (>1.4)",
                    xref="paper", yref="paper",
                    x=0.5, y=-0.15,
                    showarrow=False,
                    font=dict(size=10),
                    xanchor='center'
                )
            ]
        )
        
        # Save or show
        if output:
            self._save_figure(fig, output)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def plot_single_sample(self, result: ChromosomeAnalysisResult,
                          output: Optional[str] = None,
                          show_plot: bool = False) -> go.Figure:
        """Create bar plot showing normalized depth ratios (chromosome depth / autosomal mean).
        
        This visualization shows chromosomal copy number using normalized depth ratios,
        which is the clinical standard for CNV and aneuploidy detection. Previously
        showed absolute coverage values, which were less interpretable.
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            
        Returns:
            Plotly Figure object
        """
        # Prepare data
        chromosomes = []
        ratios = []
        colors = []
        hover_texts = []
        
        for chrom in sorted(result.chromosomes.keys(), 
                          key=ChromosomeAnalysisResult._sort_chromosomes):
            cov = result.chromosomes[chrom]
            chromosomes.append(chrom)
            ratios.append(cov.normalized_ratio)
            
            # Use biological interpretation for colors instead of status
            colors.append(self.get_ratio_color(cov.normalized_ratio))
            
            # Get biological interpretation
            interpretation = self.get_ratio_interpretation(
                cov.normalized_ratio,
                chrom,
                result.inferred_karyotype.value
            )
            
            # Enhanced hover text with actual values and interpretation
            hover_text = (
                f"<b>{chrom}</b><br>"
                f"Normalized Ratio: {cov.normalized_ratio:.3f}<br>"
                f"Mean Depth: {cov.mean_coverage:.2f}x<br>"
                f"Autosomal Mean: {result.autosomal_mean:.2f}x<br>"
                f"<b>Interpretation: {interpretation}</b><br>"
                f"Covered Bases: {cov.covered_bases:,} / {cov.total_bases:,}<br>"
                f"Percent Covered: {cov.percent_covered:.2f}%<br>"
                f"Status: {cov.status.value}"
            )
            hover_texts.append(hover_text)
        
        # Create figure
        fig = go.Figure()
        
        # Add bar chart
        fig.add_trace(go.Bar(
            x=chromosomes,
            y=ratios,
            marker_color=colors,
            hovertemplate='%{hovertext}<extra></extra>',
            hovertext=hover_texts,
            showlegend=False
        ))
        
        # Add reference lines with biological meaning
        fig.add_hline(
            y=1.0,
            line_dash="solid",
            line_color="green",
            line_width=2,
            annotation_text="Diploid (1.0)",
            annotation_position="right"
        )
        
        fig.add_hline(
            y=0.5,
            line_dash="dash",
            line_color="orange",
            line_width=1.5,
            annotation_text="Haploid (0.5)",
            annotation_position="right"
        )
        
        fig.add_hline(
            y=1.5,
            line_dash="dot",
            line_color="red",
            line_width=1.5,
            annotation_text="Trisomy (1.5)",
            annotation_position="right"
        )
        
        # Update layout
        title_text = f"Chromosome Coverage Analysis: {result.sample_id}"
        subtitle_parts = [f"Karyotype: {result.inferred_karyotype.value}"]
        subtitle_parts.append(f"Autosomal Mean: {result.autosomal_mean:.2f}x")
        
        if result.deletions:
            subtitle_parts.append(f"⚠️ Deletions: {', '.join(result.deletions)}")
        
        subtitle = " | ".join(subtitle_parts)
        
        fig.update_layout(
            title=dict(
                text=f"{title_text}<br><sub>{subtitle}</sub>",
                x=0.5,
                xanchor='center'
            ),
            xaxis_title="Chromosome",
            yaxis_title="Normalized Coverage Ratio (Chromosome Depth / Autosomal Mean)",
            width=self.width,
            height=self.height,
            template="plotly_white",
            hovermode='x',
            font=dict(size=12),
            yaxis=dict(
                range=[0, max(2.0, max(ratios) * 1.1)],
                gridcolor='#f0f0f0'
            ),
            annotations=[
                dict(
                    text="<b>Ratio Guide:</b> "
                         "<span style='color:#2ecc71'>■</span> Normal (0.9-1.1)  "
                         "<span style='color:#f39c12'>■</span> Haploid (0.4-0.6)  "
                         "<span style='color:#e67e22'>■</span> Partial Loss/Gain  "
                         "<span style='color:#e74c3c'>■</span> Deletion/Duplication",
                    xref="paper", yref="paper",
                    x=0.5, y=-0.15,
                    showarrow=False,
                    font=dict(size=10),
                    xanchor='center'
                )
            ]
        )
        
        # Save or show
        if output:
            self._save_figure(fig, output)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def plot_comparison(self, results: List[ChromosomeAnalysisResult],
                       output: Optional[str] = None,
                       show_plot: bool = False,
                       metric: str = 'ratio') -> go.Figure:
        """Create comparison plot for multiple samples.
        
        Args:
            results: List of ChromosomeAnalysisResult objects
            output: Optional output file path
            show_plot: Whether to display plot interactively
            metric: Metric to plot ('ratio' or 'coverage')
            
        Returns:
            Plotly Figure object
        """
        if not results:
            raise ValueError("No results provided")
        
        # Get all chromosomes
        all_chroms = set()
        for result in results:
            all_chroms.update(result.chromosomes.keys())
        
        chromosomes = sorted(all_chroms, key=ChromosomeAnalysisResult._sort_chromosomes)
        
        # Create figure
        fig = go.Figure()
        
        # Add trace for each sample
        for result in results:
            values = []
            hover_texts = []
            
            for chrom in chromosomes:
                if chrom in result.chromosomes:
                    cov = result.chromosomes[chrom]
                    if metric == 'ratio':
                        values.append(cov.normalized_ratio)
                        hover_text = (
                            f"<b>{result.sample_id} - {chrom}</b><br>"
                            f"Ratio: {cov.normalized_ratio:.3f}<br>"
                            f"Coverage: {cov.mean_coverage:.2f}x<br>"
                            f"Status: {cov.status.value}"
                        )
                    else:  # coverage
                        values.append(cov.mean_coverage)
                        hover_text = (
                            f"<b>{result.sample_id} - {chrom}</b><br>"
                            f"Coverage: {cov.mean_coverage:.2f}x<br>"
                            f"Ratio: {cov.normalized_ratio:.3f}<br>"
                            f"Status: {cov.status.value}"
                        )
                else:
                    values.append(0)
                    hover_text = f"<b>{result.sample_id} - {chrom}</b><br>No data"
                
                hover_texts.append(hover_text)
            
            fig.add_trace(go.Bar(
                name=result.sample_id,
                x=chromosomes,
                y=values,
                hovertemplate='%{hovertext}<extra></extra>',
                hovertext=hover_texts
            ))
        
        # Add reference lines
        if metric == 'ratio':
            fig.add_hline(y=1.0, line_dash="dash", line_color="gray",
                         annotation_text="Expected ratio: 1.0")
            fig.add_hline(y=0.15, line_dash="dot", line_color="red",
                         annotation_text="Deletion threshold")
            ylabel = "Normalized Coverage Ratio"
        else:
            ylabel = "Mean Coverage (x)"
        
        # Update layout
        fig.update_layout(
            title=f"Multi-Sample Chromosome Coverage Comparison ({len(results)} samples)",
            xaxis_title="Chromosome",
            yaxis_title=ylabel,
            width=self.width,
            height=self.height,
            template="plotly_white",
            barmode='group',
            hovermode='x unified',
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            )
        )
        
        # Save or show
        if output:
            self._save_figure(fig, output)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def plot_heatmap(self, results: List[ChromosomeAnalysisResult],
                    output: Optional[str] = None,
                    show_plot: bool = False) -> go.Figure:
        """Create heatmap of normalized coverage ratios across samples.
        
        Args:
            results: List of ChromosomeAnalysisResult objects
            output: Optional output file path
            show_plot: Whether to display plot interactively
            
        Returns:
            Plotly Figure object
        """
        if not results:
            raise ValueError("No results provided")
        
        # Get all chromosomes
        all_chroms = set()
        for result in results:
            all_chroms.update(result.chromosomes.keys())
        
        chromosomes = sorted(all_chroms, key=ChromosomeAnalysisResult._sort_chromosomes)
        sample_ids = [r.sample_id for r in results]
        
        # Build matrix
        matrix = []
        hover_texts = []
        
        for result in results:
            row = []
            hover_row = []
            for chrom in chromosomes:
                if chrom in result.chromosomes:
                    cov = result.chromosomes[chrom]
                    row.append(cov.normalized_ratio)
                    hover_text = (
                        f"Sample: {result.sample_id}<br>"
                        f"Chromosome: {chrom}<br>"
                        f"Ratio: {cov.normalized_ratio:.3f}<br>"
                        f"Coverage: {cov.mean_coverage:.2f}x<br>"
                        f"Status: {cov.status.value}"
                    )
                else:
                    row.append(0)
                    hover_text = f"Sample: {result.sample_id}<br>Chromosome: {chrom}<br>No data"
                hover_row.append(hover_text)
            matrix.append(row)
            hover_texts.append(hover_row)
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=chromosomes,
            y=sample_ids,
            colorscale='RdYlGn',
            zmid=1.0,  # Center at expected ratio
            zmin=0,
            zmax=2,
            hovertemplate='%{customdata}<extra></extra>',
            customdata=hover_texts,
            colorbar=dict(title="Normalized<br>Ratio")
        ))
        
        fig.update_layout(
            title="Chromosome Coverage Heatmap",
            xaxis_title="Chromosome",
            yaxis_title="Sample",
            width=self.width,
            height=max(400, len(results) * 40),
            template="plotly_white"
        )
        
        # Save or show
        if output:
            self._save_figure(fig, output)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def plot_clinical_report(self, result: ChromosomeAnalysisResult,
                            output: Optional[str] = None,
                            show_plot: bool = False) -> go.Figure:
        """Create two-panel clinical report for WES data - the recommended comprehensive view.
        
        This creates a publication-quality clinical report with two vertically stacked panels:
        
        **Top Panel (Clinical Signal)**: Normalized depth ratio
        - Primary metric for CNV/aneuploidy detection
        - Ratio = chromosome depth / autosomal mean
        - Reference lines at 1.0 (diploid), 0.5 (haploid), 1.5 (trisomy)
        - Colors based on biological interpretation
        
        **Bottom Panel (Data Quality)**: Absolute coverage depth
        - Shows actual sequencing depth per chromosome
        - Helps assess data quality and coverage uniformity
        - Reference line at autosomal mean
        
        Both panels share the same X-axis (chromosomes) for easy comparison.
        Optimized for WES data interpretation.
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            
        Returns:
            Plotly Figure object
        """
        # Prepare data
        chromosomes = sorted(result.chromosomes.keys(),
                           key=ChromosomeAnalysisResult._sort_chromosomes)
        
        ratios = []
        coverages = []
        colors_ratio = []
        colors_coverage = []
        hover_texts_ratio = []
        hover_texts_coverage = []
        
        for chrom in chromosomes:
            cov = result.chromosomes[chrom]
            ratios.append(cov.normalized_ratio)
            coverages.append(cov.mean_coverage)
            
            # Colors based on biological interpretation for ratio
            colors_ratio.append(self.get_ratio_color(cov.normalized_ratio))
            
            # Colors for coverage based on quality (relative to autosomal mean)
            if cov.mean_coverage >= result.autosomal_mean * 0.8:
                colors_coverage.append('#3498db')  # Blue - good depth
            elif cov.mean_coverage >= result.autosomal_mean * 0.5:
                colors_coverage.append('#f39c12')  # Orange - moderate depth
            else:
                colors_coverage.append('#e74c3c')  # Red - low depth
            
            # Get interpretation
            interpretation = self.get_ratio_interpretation(
                cov.normalized_ratio,
                chrom,
                result.inferred_karyotype.value
            )
            
            # Hover text for ratio panel
            hover_texts_ratio.append(
                f"<b>{chrom}</b><br>"
                f"Normalized Ratio: {cov.normalized_ratio:.3f}<br>"
                f"<b>Interpretation: {interpretation}</b><br>"
                f"Mean Depth: {cov.mean_coverage:.2f}x<br>"
                f"Autosomal Mean: {result.autosomal_mean:.2f}x"
            )
            
            # Hover text for coverage panel
            hover_texts_coverage.append(
                f"<b>{chrom}</b><br>"
                f"Mean Depth: {cov.mean_coverage:.2f}x<br>"
                f"Autosomal Mean: {result.autosomal_mean:.2f}x<br>"
                f"Covered Bases: {cov.covered_bases:,}<br>"
                f"Percent Covered: {cov.percent_covered:.2f}%"
            )
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=(
                '<b>Clinical Signal:</b> Normalized Depth Ratio (CNV/Aneuploidy Detection)',
                '<b>Data Quality:</b> Absolute Coverage Depth'
            ),
            vertical_spacing=0.12,
            row_heights=[0.55, 0.45]
        )
        
        # Top panel: Normalized ratios (clinical signal)
        fig.add_trace(
            go.Bar(
                x=chromosomes,
                y=ratios,
                marker_color=colors_ratio,
                name="Normalized Ratio",
                hovertemplate='%{hovertext}<extra></extra>',
                hovertext=hover_texts_ratio,
                showlegend=False
            ),
            row=1, col=1
        )
        
        # Add reference lines to ratio panel
        fig.add_hline(
            y=1.0, line_dash="solid", line_color="green", line_width=2,
            annotation_text="Diploid (1.0)", annotation_position="right",
            row=1, col=1
        )
        fig.add_hline(
            y=0.5, line_dash="dash", line_color="orange", line_width=1.5,
            annotation_text="Haploid (0.5)", annotation_position="right",
            row=1, col=1
        )
        fig.add_hline(
            y=1.5, line_dash="dot", line_color="red", line_width=1.5,
            annotation_text="Trisomy (1.5)", annotation_position="right",
            row=1, col=1
        )
        
        # Bottom panel: Absolute coverage (data quality)
        fig.add_trace(
            go.Bar(
                x=chromosomes,
                y=coverages,
                marker_color=colors_coverage,
                name="Coverage Depth",
                hovertemplate='%{hovertext}<extra></extra>',
                hovertext=hover_texts_coverage,
                showlegend=False
            ),
            row=2, col=1
        )
        
        # Add reference line at autosomal mean
        fig.add_hline(
            y=result.autosomal_mean, line_dash="dash", line_color="gray",
            annotation_text=f"Autosomal Mean: {result.autosomal_mean:.1f}x",
            annotation_position="right",
            row=2, col=1
        )
        
        # Update axes
        fig.update_xaxes(title_text="Chromosome", row=2, col=1)
        fig.update_yaxes(
            title_text="Normalized Ratio",
            range=[0, max(2.0, max(ratios) * 1.1)],
            gridcolor='#f0f0f0',
            row=1, col=1
        )
        fig.update_yaxes(
            title_text="Mean Coverage (x)",
            gridcolor='#f0f0f0',
            row=2, col=1
        )
        
        # Update layout
        title_parts = [f"Sample: {result.sample_id}"]
        title_parts.append(f"Karyotype: {result.inferred_karyotype.value}")
        title_parts.append(f"Autosomal Mean: {result.autosomal_mean:.2f}x")
        
        if result.deletions:
            title_parts.append(f"⚠️ Deletions: {', '.join(result.deletions)}")
        
        subtitle = " | ".join(title_parts)
        
        fig.update_layout(
            title=dict(
                text=f"<b>Clinical WES Report: Chromosome Coverage Analysis</b><br>"
                     f"<sub>{subtitle}</sub>",
                x=0.5,
                xanchor='center',
                font=dict(size=16)
            ),
            width=self.width,
            height=self.height * 1.5,
            template="plotly_white",
            hovermode='x',
            font=dict(size=11),
            showlegend=False,
            annotations=[
                dict(
                    text="<b>Top panel:</b> Use ratios to detect CNVs/aneuploidies | "
                         "<b>Bottom panel:</b> Verify data quality and coverage depth",
                    xref="paper", yref="paper",
                    x=0.5, y=-0.08,
                    showarrow=False,
                    font=dict(size=10),
                    xanchor='center'
                )
            ]
        )
        
        # Save or show
        if output:
            self._save_figure(fig, output)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def create_report_plot(self, result: ChromosomeAnalysisResult,
                          output: Optional[str] = None) -> go.Figure:
        """Create comprehensive report with multiple subplots.
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path
            
        Returns:
            Plotly Figure object
        """
        # Create subplots
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=(
                'Absolute Coverage per Chromosome',
                'Normalized Coverage Ratios'
            ),
            vertical_spacing=0.15,
            row_heights=[0.5, 0.5]
        )
        
        # Prepare data
        chromosomes = sorted(result.chromosomes.keys(),
                           key=ChromosomeAnalysisResult._sort_chromosomes)
        
        coverages = []
        ratios = []
        colors = []
        
        for chrom in chromosomes:
            cov = result.chromosomes[chrom]
            coverages.append(cov.mean_coverage)
            ratios.append(cov.normalized_ratio)
            colors.append(self.COLORS[cov.status])
        
        # Plot 1: Absolute coverage
        fig.add_trace(
            go.Bar(x=chromosomes, y=coverages, marker_color=colors,
                  name="Coverage", showlegend=False),
            row=1, col=1
        )
        fig.add_hline(y=result.autosomal_mean, line_dash="dash",
                     line_color="gray", row=1, col=1)
        
        # Plot 2: Normalized ratios
        fig.add_trace(
            go.Bar(x=chromosomes, y=ratios, marker_color=colors,
                  name="Ratio", showlegend=False),
            row=2, col=1
        )
        fig.add_hline(y=1.0, line_dash="dash", line_color="gray", row=2, col=1)
        fig.add_hline(y=0.15, line_dash="dot", line_color="red", row=2, col=1)
        
        # Update axes
        fig.update_xaxes(title_text="Chromosome", row=2, col=1)
        fig.update_yaxes(title_text="Coverage (x)", row=1, col=1)
        fig.update_yaxes(title_text="Normalized Ratio", row=2, col=1)
        
        # Update layout
        title_text = f"Chromosome Analysis Report: {result.sample_id}"
        if result.inferred_karyotype:
            title_text += f" | Karyotype: {result.inferred_karyotype.value}"
        
        fig.update_layout(
            title=title_text,
            width=self.width,
            height=self.height * 1.5,
            template="plotly_white",
            showlegend=False
        )
        
        if output:
            self._save_figure(fig, output)
        
        return fig
    
    def _save_figure(self, fig: go.Figure, output: str):
        """Save figure to file.
        
        Args:
            fig: Plotly figure
            output: Output file path
        """
        output_path = Path(output)
        
        if output_path.suffix.lower() == '.html':
            fig.write_html(str(output_path))
        elif output_path.suffix.lower() in ['.png', '.jpg', '.jpeg', '.svg', '.pdf']:
            # Requires kaleido package
            try:
                fig.write_image(str(output_path), width=self.width, height=self.height)
            except Exception as e:
                # Fallback to HTML if image export fails
                html_path = output_path.with_suffix('.html')
                fig.write_html(str(html_path))
                print(f"Warning: Could not save as {output_path.suffix}. "
                      f"Saved as HTML instead: {html_path}")
                print(f"Install kaleido for image export: pip install kaleido")
        else:
            # Default to HTML
            fig.write_html(str(output_path))
    
    @staticmethod
    def generate_ascii_bar(value: float, max_value: float, width: int = 40) -> str:
        """Generate ASCII bar for text reports.
        
        Args:
            value: Current value
            max_value: Maximum value for scaling
            width: Bar width in characters
            
        Returns:
            ASCII bar string
        """
        if max_value <= 0:
            return "▌" + " " * (width - 1)
        
        filled = int((value / max_value) * width)
        filled = max(0, min(filled, width))
        
        bar = "█" * filled + "░" * (width - filled)
        return bar
