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
    
    def plot_chromosome_tracks(self, result: ChromosomeAnalysisResult,
                              output: Optional[str] = None,
                              show_plot: bool = False,
                              transformation: str = 'log') -> go.Figure:
        """Create simplified chromosome coverage tracks with percentage-based fills.
        
        This visualization shows each chromosome as a horizontal bar where:
        - Gray background = full chromosome (100%)
        - Colored fill = percentage of chromosome covered (log-scaled for visibility)
        - Vertical line = expected WES coverage for that chromosome
        
        Uses logarithmic transformation to make small coverages (even single genes) visible.
        Perfect for WES/targeted sequencing to detect deletions and aneuploidies.
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            transformation: Transform method ('log', 'sqrt', 'linear') for visibility
            
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
            
            # Get actual and transformed coverage percentages
            actual_percent = cov.percent_covered
            visual_percent = self._transform_coverage_for_visibility(actual_percent, transformation)
            
            # Get color based on status
            color = self.COLORS.get(cov.status, '#95a5a6')
            
            # Add gray background (full chromosome = 100%)
            all_shapes.append(dict(
                type="rect",
                x0=0, x1=100,
                y0=y_pos - 0.4, y1=y_pos + 0.4,
                fillcolor='#ecf0f1',
                line=dict(color='#bdc3c7', width=1),
                layer='below'
            ))
            
            # Add colored fill for covered percentage (transformed for visibility)
            if visual_percent > 0:
                all_shapes.append(dict(
                    type="rect",
                    x0=0, x1=visual_percent,
                    y0=y_pos - 0.4, y1=y_pos + 0.4,
                    fillcolor=color,
                    line=dict(color=color, width=0),
                    opacity=0.8
                ))
            
            # Add expected WES coverage line for this chromosome
            expected_cov = self.EXPECTED_WES_COVERAGE.get(chrom, 2.0)
            expected_visual = self._transform_coverage_for_visibility(expected_cov, transformation)
            all_shapes.append(dict(
                type="line",
                x0=expected_visual, x1=expected_visual,
                y0=y_pos - 0.5, y1=y_pos + 0.5,
                line=dict(color='rgba(0,0,0,0.3)', width=1, dash='dot')
            ))
            
            # Add invisible scatter trace for hover information (one point per chromosome)
            hover_text = (
                f"<b>{chrom}</b><br>"
                f"Actual coverage: {actual_percent:.3f}% of chromosome<br>"
                f"Visual (log-scaled): {visual_percent:.1f}%<br>"
                f"Expected WES: ~{expected_cov:.1f}%<br>"
                f"Mean depth: {cov.mean_coverage:.2f}x<br>"
                f"Status: {cov.status.value}<br>"
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
        
        # Count chromosomes by status for summary
        status_counts = {}
        for cov in result.chromosomes.values():
            status = cov.status.value
            status_counts[status] = status_counts.get(status, 0) + 1
        
        # Create status summary for subtitle
        status_summary = " | ".join([f"{count} {status}" for status, count in sorted(status_counts.items())])
        
        # Transformation note for subtitle
        transform_note = {
            'log': 'Log-scaled for visibility',
            'sqrt': 'Square-root scaled for visibility',
            'linear': 'Linear scale'
        }.get(transformation, 'Transformed for visibility')
        
        # Update layout with proper margins and axis settings
        fig.update_layout(
            title=dict(
                text=f"Chromosome Coverage Map: {result.sample_id}<br>"
                     f"<sub>Karyotype: {result.inferred_karyotype.value} | "
                     f"Autosomal mean: {result.autosomal_mean:.2f}x<br>"
                     f"{status_summary}<br>"
                     f"<i>{transform_note} | Dotted lines show expected WES coverage | Hover for actual percentages</i></sub>",
                x=0.5,
                xanchor='center'
            ),
            xaxis=dict(
                title="Coverage (% of chromosome, transformed for visibility)",
                showgrid=True,
                gridcolor='#f0f0f0',
                zeroline=False,
                range=[0, 105]  # 0-100% with small padding
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
            margin=dict(l=80, r=40, t=140, b=80)  # Left margin for chromosome labels
        )
        
        # Save if output specified
        if output:
            output_path = Path(output)
            if output_path.suffix == '.html':
                fig.write_html(str(output_path))
            else:
                fig.write_image(str(output_path), width=self.width, 
                              height=self.height, scale=self.dpi/100)
        
        if show_plot:
            fig.show()
        
        return fig
    
    def plot_single_sample(self, result: ChromosomeAnalysisResult,
                          output: Optional[str] = None,
                          show_plot: bool = False) -> go.Figure:
        """Create bar plot for a single sample's chromosome coverage.
        
        Args:
            result: ChromosomeAnalysisResult to visualize
            output: Optional output file path (HTML or image)
            show_plot: Whether to display plot interactively
            
        Returns:
            Plotly Figure object
        """
        # Prepare data
        chromosomes = []
        coverages = []
        ratios = []
        statuses = []
        colors = []
        hover_texts = []
        
        for chrom in sorted(result.chromosomes.keys(), 
                          key=ChromosomeAnalysisResult._sort_chromosomes):
            cov = result.chromosomes[chrom]
            chromosomes.append(chrom)
            coverages.append(cov.mean_coverage)
            ratios.append(cov.normalized_ratio)
            
            # Create status and hover text based on chromosome state
            statuses.append(cov.status.value)
            colors.append(self.COLORS[cov.status])
            
            if cov.status == ChromosomeStatus.NOT_SEQUENCED:
                hover_text = (
                    f"<b>{chrom}</b><br>"
                    f"<i>Not sequenced - no data in BAM</i><br>"
                    f"This chromosome was not included in the sequencing"
                )
            elif cov.status == ChromosomeStatus.PARTIAL:
                hover_text = (
                    f"<b>{chrom}</b><br>"
                    f"<i>Partially sequenced ({cov.percent_covered:.1f}% of chromosome)</i><br>"
                    f"Mean Coverage: {cov.mean_coverage:.2f}x (in sequenced regions)<br>"
                    f"Normalized Ratio: {cov.normalized_ratio:.3f}<br>"
                    f"Status: {cov.status.value}<br>"
                    f"Note: Only targeted regions sequenced (e.g., WES)"
                )
            else:
                hover_text = (
                    f"<b>{chrom}</b><br>"
                    f"Mean Coverage: {cov.mean_coverage:.2f}x<br>"
                    f"Normalized Ratio: {cov.normalized_ratio:.3f}<br>"
                    f"Status: {cov.status.value}<br>"
                    f"Covered: {cov.covered_bases:,} / {cov.total_bases:,} bp<br>"
                    f"Chromosome coverage: {cov.percent_covered:.1f}%"
                )
            
            hover_texts.append(hover_text)
        
        # Create figure
        fig = go.Figure()
        
        # Add bar chart
        fig.add_trace(go.Bar(
            x=chromosomes,
            y=coverages,
            marker_color=colors,
            hovertemplate='%{hovertext}<extra></extra>',
            hovertext=hover_texts,
            showlegend=False
        ))
        
        # Add reference line at autosomal mean
        fig.add_hline(
            y=result.autosomal_mean,
            line_dash="dash",
            line_color="gray",
            annotation_text=f"Autosomal Mean: {result.autosomal_mean:.1f}x",
            annotation_position="right"
        )
        
        # Update layout
        title_text = f"Chromosome Coverage Analysis: {result.sample_id}"
        if result.inferred_karyotype:
            title_text += f"<br><sub>Inferred Karyotype: {result.inferred_karyotype.value}</sub>"
        
        # Count chromosomes by status for subtitle
        sequenced_count = sum(1 for c in result.chromosomes.values() if c.mean_coverage > 0)
        total_count = len(result.chromosomes)
        
        fig.update_layout(
            title=title_text,
            xaxis_title="Chromosome",
            yaxis_title="Mean Coverage (x)",
            width=self.width,
            height=self.height,
            template="plotly_white",
            hovermode='x',
            font=dict(size=12),
            annotations=[
                dict(
                    text=f"Chromosomes sequenced: {sequenced_count}/{total_count} | "
                         f"<span style='color:#2ecc71'>■</span> Normal  "
                         f"<span style='color:#9b59b6'>■</span> Partial (WES)  "
                         f"<span style='color:#f39c12'>■</span> Reduced  "
                         f"<span style='color:#e74c3c'>■</span> Deleted  "
                         f"<span style='color:#95a5a6'>■</span> Not Sequenced",
                    xref="paper", yref="paper",
                    x=0.5, y=-0.15,
                    showarrow=False,
                    font=dict(size=10),
                    xanchor='center'
                )
            ]
        )
        
        # Add warnings as annotations if present
        if result.deletions:
            warning_text = "⚠️ Deletions detected: " + ", ".join(result.deletions)
            fig.add_annotation(
                text=warning_text,
                xref="paper", yref="paper",
                x=0.5, y=1.08,
                showarrow=False,
                font=dict(color="red", size=14),
                xanchor='center'
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
