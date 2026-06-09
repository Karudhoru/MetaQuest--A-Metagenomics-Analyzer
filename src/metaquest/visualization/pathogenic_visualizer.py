#!/usr/bin/env python3
"""
MetaQuest Visualization Module
Pathogenic Visualizer Class - Fixed for integrated risk reports
"""
from collections import Counter
import json
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .base_visualizer import BaseVisualizer

class PathogenVisualizer(BaseVisualizer):
    """Pathogen detection visualizations aligned with reporting"""
    
    def _create_summary_risk_chart(self, report_data):
        """Create a summary chart when only aggregate data is available"""
        try:
            self.formatter.info("Creating summary risk chart from aggregate data")
            
            counts = report_data.get('pathogen_counts', {})
            markers = report_data.get('functional_markers', {})
            scores = report_data.get('risk_scores', {})
            
            # Create a bar chart showing risk components
            categories = []
            values = []
            colors = []
            hover_text = []
            
            # Add pathogen counts
            if counts.get('high_risk', 0) > 0:
                categories.append("High Risk Pathogens")
                values.append(counts['high_risk'])
                colors.append('#DC143C')
                hover_text.append(f"High Risk: {counts['high_risk']} pathogens detected")
            
            if counts.get('moderate_risk', 0) > 0:
                categories.append("Moderate Risk Pathogens")
                values.append(counts['moderate_risk'])
                colors.append('#FFA500')
                hover_text.append(f"Moderate Risk: {counts['moderate_risk']} pathogens detected")
            
            if counts.get('low_risk', 0) > 0:
                categories.append("Low Risk Pathogens")
                values.append(counts['low_risk'])
                colors.append('#32CD32')
                hover_text.append(f"Low Risk: {counts['low_risk']} pathogens detected")
            
            # Add functional markers (scaled for visualization)
            if markers.get('virulence', 0) > 0:
                categories.append("Virulence Genes")
                scaled_val = min(markers['virulence'] / 20, 10)
                values.append(scaled_val)
                colors.append('#8B0000')
                hover_text.append(f"Virulence: {markers['virulence']} genes detected")
            
            if markers.get('amr', 0) > 0:
                categories.append("AMR Genes")
                scaled_val = min(markers['amr'] / 30, 10)
                values.append(scaled_val)
                colors.append('#FF4500')
                hover_text.append(f"AMR: {markers['amr']} genes detected")
            
            if markers.get('transposases', 0) > 0:
                categories.append("Transposases")
                scaled_val = min(markers['transposases'] / 5, 10)
                values.append(scaled_val)
                colors.append('#FFD700')
                hover_text.append(f"Transposases: {markers['transposases']} genes detected")
            
            if not categories:
                self.formatter.warning("No data available for summary chart") 
                return None
            
            fig = go.Figure(data=[go.Bar(
                x=categories,
                y=values,
                marker_color=colors,
                text=[f"{v:.1f}" for v in values],
                textposition='outside',
                hovertext=hover_text,
                hoverinfo='text'
            )])
            
            risk_level = scores.get('risk_level', 'Unknown')
            final_score = scores.get('final', 0)
            total_pathogens = counts.get('total', 0)
            
            fig.update_layout(
                title=f"Integrated Pathogen Risk Summary<br><sub>Overall Risk: <b>{risk_level}</b> (Score: {final_score:.1f}/100) | {total_pathogens} Pathogens Detected</sub>",
                xaxis_title="Risk Component",
                yaxis_title="Count / Scaled Value",
                template="plotly_white",
                height=500,
                xaxis={'tickangle': -45},
                margin=dict(b=150, t=100),
                font=dict(size=11),
                showlegend=False
            )
            
            filepath = self.save_plot(fig, "pathogen_risk_summary.html")
            self.formatter.success(f"Risk summary chart saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating summary chart: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_risk_assessment_chart(self, report_path, title="Pathogen Risk Assessment"):
        """
        Create comprehensive risk assessment visualization.
        Supports FASTQ, FASTA+ML, and integrated summary formats.
        """
        try:
            # Load report data
            report_data = self.load_json_report(Path(report_path))
            if not report_data:
                return None
            
            # Check if this is an integrated risk summary (no detailed detections)
            if 'pathogen_counts' in report_data and 'data' not in report_data:
                self.formatter.info("Detected integrated risk summary format")
                return self._create_summary_risk_chart(report_data)
            
            # Extract the nested data structure
            data_section = report_data.get('data', {})
            summary = data_section.get('summary', {})
            
            # Determine report type and extract detections
            analysis_type = summary.get('analysis_type', 'unknown')
            
            # Try multiple paths to find pathogen data
            if 'FASTA' in analysis_type or 'BLAST' in analysis_type:
                detections = data_section.get('blast_taxonomy', {}).get('organisms', {})
                self.formatter.info("Processing FASTA+ML pathogen report")
            else:
                # For FASTQ reports, try multiple possible locations
                detections = (data_section.get('pathogen_detections') or 
                            data_section.get('pathogens') or
                            report_data.get('pathogen_detections') or
                            report_data.get('pathogens') or {})
                self.formatter.info("Processing FASTQ pathogen report")
            
            if not detections:
                self.formatter.warning("No pathogen detections found in report")
                self.formatter.debug(f"Available keys in report root: {list(report_data.keys())}")
                self.formatter.debug(f"Available keys in data section: {list(data_section.keys())}")
                return None

            # Convert to DataFrame
            data_rows = []
            for organism, details in detections.items():
                # Handle different data structures
                risk_level = details.get('risk_level', 'Unknown')
                abundance = details.get('abundance_percentage', 0)
                
                # Get identity from various possible keys
                identity = details.get('sequence_identity') or \
                          details.get('max_identity') or \
                          details.get('blast_identity', 0)
                
                hits = details.get('sequence_hits', 0) or \
                      details.get('hits', 0) or \
                      details.get('blast_hits', 0)
                
                confidence = details.get('confidence_score', 0)
                methods = details.get('detection_methods', [])
                method_count = len(methods) if isinstance(methods, list) else 0
                
                data_rows.append({
                    'organism': organism,
                    'risk_level': risk_level,
                    'abundance': abundance,
                    'identity': identity,
                    'hits': hits,
                    'confidence': confidence,
                    'method_count': method_count
                })

            df = pd.DataFrame(data_rows)
            if df.empty:
                self.formatter.warning("No pathogen data to visualize")
                return None

            self.formatter.info(f"Loaded {len(df)} pathogens for risk assessment chart")

            # Sort by risk and relevant metric
            risk_order = {'CRITICAL': 4, 'HIGH': 3, 'MEDIUM': 2, 'LOW': 1, 'Unknown': 0}
            df['risk_numeric'] = df['risk_level'].map(risk_order)
            df = df.sort_values(['risk_numeric', 'hits'], ascending=[False, False])
            top_pathogens = df.head(20)

            # Create visualization
            risk_colors = {
                'CRITICAL': '#DC143C',
                'HIGH': '#FF6347',
                'MEDIUM': '#FFA500',
                'LOW': '#32CD32',
                'Unknown': '#808080'
            }
            
            # Use appropriate metric for y-axis
            if top_pathogens['hits'].sum() > 0:
                y_values = top_pathogens['hits']
                y_label = "Detection Count"
            elif top_pathogens['abundance'].sum() > 0:
                y_values = top_pathogens['abundance']
                y_label = "Abundance (%)"
            else:
                y_values = pd.Series([1] * len(top_pathogens))
                y_label = "Presence"

            fig = go.Figure()
            
            for risk_level in ['CRITICAL', 'HIGH', 'MEDIUM', 'LOW', 'Unknown']:
                mask = top_pathogens['risk_level'] == risk_level
                if mask.any():
                    subset = top_pathogens[mask]
                    fig.add_trace(go.Bar(
                        name=risk_level,
                        x=subset['organism'],
                        y=y_values[mask],
                        marker_color=risk_colors[risk_level],
                        hovertemplate=(
                            '<b>%{x}</b><br>' +
                            f'Risk: {risk_level}<br>' +
                            'Count: %{y}<br>' +
                            'Identity: %{customdata[0]:.1f}%<br>' +
                            'Confidence: %{customdata[1]:.3f}<br>' +
                            'Methods: %{customdata[2]}<extra></extra>'
                        ),
                        customdata=list(zip(subset['identity'], subset['confidence'], subset['method_count']))
                    ))
            
            overall_risk = summary.get('overall_risk_assessment', 'Unknown')
            fig.update_layout(
                title=f"{title}<br><sub>Overall Risk Level: <b>{overall_risk}</b></sub>",
                xaxis_title="Pathogen Species",
                yaxis_title=y_label,
                template="plotly_white",
                height=600,
                xaxis={'tickangle': -45},
                margin=dict(b=180, t=100),
                font=dict(size=11),
                legend=dict(
                    title="Risk Level",
                    orientation="v",
                    yanchor="top",
                    y=0.99,
                    xanchor="right",
                    x=0.99
                ),
                barmode='group'
            )

            filepath = self.save_plot(fig, "pathogen_risk_assessment.html")
            self.formatter.success(f"Risk assessment chart saved: {Path(filepath).name}")
            return filepath

        except Exception as e:
            self.formatter.error(f"Error creating risk assessment chart: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_who_priority_distribution(self, report_path):
        """Create WHO priority pathogen distribution chart"""
        try:
            report_data = self.load_json_report(Path(report_path))
            if not report_data:
                return None
            
            # Check for integrated summary format
            if 'pathogen_counts' in report_data and 'data' not in report_data:
                self.formatter.info("Skipping WHO chart - integrated summary lacks WHO priority data")
                return None
            
            # Extract nested data
            data_section = report_data.get('data', {})
            summary = data_section.get('summary', {})
            analysis_type = summary.get('analysis_type', '')
            
            # Extract detections based on analysis type
            if 'FASTA' in analysis_type or 'BLAST' in analysis_type:
                detections = data_section.get('blast_taxonomy', {}).get('organisms', {})
            else:
                detections = (data_section.get('pathogen_detections') or 
                            data_section.get('pathogens') or
                            report_data.get('pathogen_detections') or
                            report_data.get('pathogens') or {})
            
            if not detections:
                self.formatter.warning("No detections for WHO priority chart")
                return None
            
            # Count by WHO priority
            who_counts = Counter()
            for organism, details in detections.items():
                priority = details.get('who_priority', 'Not Listed')
                who_counts[priority] += 1
            
            if not who_counts:
                self.formatter.warning("No WHO priority data to plot")
                return None
            
            self.formatter.info(f"WHO priority breakdown: {dict(who_counts)}")
            
            # Create pie chart
            labels = list(who_counts.keys())
            values = list(who_counts.values())
            
            colors = {
                'Critical Priority': '#DC143C',
                'High Priority': '#FF6347',
                'Medium Priority': '#FFA500',
                'Not Listed': '#808080'
            }
            color_list = [colors.get(label, '#808080') for label in labels]
            
            fig = go.Figure(data=[go.Pie(
                labels=labels,
                values=values,
                hole=0.4,
                marker_colors=color_list,
                hovertemplate="<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>",
                textinfo='label+percent',
                textposition='outside'
            )])
            
            fig.update_layout(
                title="WHO Priority Pathogen Classification<br><sub>Distribution by WHO 2024 Priority List</sub>",
                template="plotly_white",
                height=500,
                font=dict(size=11),
                annotations=[dict(
                    text=f"Total<br>Pathogens<br>{sum(values)}",
                    x=0.5, y=0.5,
                    font_size=16,
                    showarrow=False
                )],
                showlegend=True
            )
            
            filepath = self.save_plot(fig, "who_priority_distribution.html")
            self.formatter.success(f"WHO priority chart saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating WHO priority chart: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_detection_confidence_chart(self, report_path):
        """Create detection confidence visualization"""
        try:
            report_data = self.load_json_report(Path(report_path))
            if not report_data:
                return None
            
            # Check for integrated summary format
            if 'pathogen_counts' in report_data and 'data' not in report_data:
                self.formatter.info("Skipping confidence chart - integrated summary lacks detailed detection data")
                return None
            
            # Extract nested data
            data_section = report_data.get('data', {})
            summary = data_section.get('summary', {})
            analysis_type = summary.get('analysis_type', '')
            
            # Extract detections
            if 'FASTA' in analysis_type:
                detections = data_section.get('blast_taxonomy', {}).get('organisms', {})
            else:
                detections = (data_section.get('pathogen_detections') or 
                            data_section.get('pathogens') or
                            report_data.get('pathogen_detections') or
                            report_data.get('pathogens') or {})
            
            if not detections:
                self.formatter.warning("No detections for confidence chart")
                return None
            
            # Prepare data for scatter plot
            organisms = []
            confidences = []
            risk_levels = []
            detection_methods = []
            abundances = []
            
            for organism, details in detections.items():
                conf = details.get('confidence_score', 0)
                # Include all organisms, even if confidence is 0
                organisms.append(organism[:40])  # Truncate long names
                confidences.append(conf if conf > 0 else 0.01)  # Minimum value for visibility
                risk_levels.append(details.get('risk_level', 'Unknown'))
                methods = details.get('detection_methods', [])
                detection_methods.append(len(methods) if isinstance(methods, list) else 1)
                abundances.append(details.get('abundance_percentage', 0))
            
            if not confidences:
                self.formatter.warning("No confidence scores available")
                return None
            
            df = pd.DataFrame({
                'organism': organisms,
                'confidence': confidences,
                'risk_level': risk_levels,
                'method_count': detection_methods,
                'abundance': abundances
            })
            
            # Sort by confidence descending
            df = df.sort_values('confidence', ascending=False).head(20)
            
            self.formatter.info(f"Creating confidence chart for {len(df)} pathogens")
            
            # Create scatter plot
            risk_colors = {
                'CRITICAL': '#DC143C',
                'HIGH': '#FF6347',
                'MEDIUM': '#FFA500',
                'LOW': '#32CD32',
                'Unknown': '#808080'
            }
            
            fig = px.scatter(
                df,
                x='confidence',
                y='organism',
                color='risk_level',
                size='method_count',
                color_discrete_map=risk_colors,
                title="Pathogen Detection Confidence Analysis<br><sub>Confidence scores by detection methods</sub>",
                labels={
                    'confidence': 'Confidence Score',
                    'organism': 'Pathogen',
                    'method_count': 'Methods',
                    'risk_level': 'Risk Level'
                },
                hover_data={
                    'confidence': ':.3f', 
                    'method_count': True,
                    'abundance': ':.4f',
                    'risk_level': True
                }
            )
            
            # Add confidence threshold lines
            fig.add_vline(x=0.7, line_dash="dash", line_color="green",
                         annotation_text="High Confidence", annotation_position="top")
            fig.add_vline(x=0.5, line_dash="dash", line_color="orange",
                         annotation_text="Medium Confidence", annotation_position="top")
            fig.add_vline(x=0.3, line_dash="dash", line_color="red",
                         annotation_text="Low Confidence", annotation_position="top")
            
            fig.update_layout(
                template="plotly_white",
                height=max(500, len(df) * 30),
                xaxis={'range': [0, 1.05]},
                margin=dict(l=250, r=100, t=100, b=60),
                font=dict(size=10),
                showlegend=True,
                legend=dict(
                    title="Risk Level",
                    orientation="v",
                    yanchor="top",
                    y=0.99,
                    xanchor="right",
                    x=0.99
                )
            )
            
            filepath = self.save_plot(fig, "detection_confidence.html")
            self.formatter.success(f"Confidence chart saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating confidence chart: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_abundance_heatmap(self, report_path):
        """Create abundance heatmap for top pathogens"""
        try:
            report_data = self.load_json_report(Path(report_path))
            if not report_data:
                return None
            
            # Check for integrated summary format
            if 'pathogen_counts' in report_data and 'data' not in report_data:
                self.formatter.debug("Skipping abundance heatmap - integrated summary format")
                return None
            
            # Extract data
            data_section = report_data.get('data', {})
            detections = (data_section.get('pathogen_detections') or 
                         data_section.get('pathogens') or
                         report_data.get('pathogen_detections') or
                         report_data.get('pathogens') or {})
            
            if not detections:
                self.formatter.debug("No detections for abundance heatmap")
                return None
            
            # Filter pathogens with abundance data
            pathogen_data = []
            for organism, details in detections.items():
                abundance = details.get('abundance_percentage', 0)
                if abundance > 0:
                    pathogen_data.append({
                        'organism': organism[:35],
                        'abundance': abundance,
                        'risk_level': details.get('risk_level', 'Unknown'),
                        'confidence': details.get('confidence_score', 0)
                    })
            
            if not pathogen_data:
                self.formatter.debug("No abundance data available for heatmap")
                return None
            
            df = pd.DataFrame(pathogen_data)
            df = df.sort_values('abundance', ascending=False).head(15)
            
            self.formatter.debug(f"Creating abundance heatmap for {len(df)} pathogens")
            
            # Create horizontal bar chart with color gradient
            fig = go.Figure()
            
            # Color scale based on risk
            risk_colors = {
                'CRITICAL': '#DC143C',
                'HIGH': '#FF6347',
                'MEDIUM': '#FFA500',
                'LOW': '#32CD32',
                'Unknown': '#808080'
            }
            
            colors = [risk_colors.get(risk, '#808080') for risk in df['risk_level']]
            
            fig.add_trace(go.Bar(
                x=df['abundance'],
                y=df['organism'],
                orientation='h',
                marker=dict(
                    color=colors,
                    line=dict(color='black', width=1)
                ),
                text=df['abundance'].apply(lambda x: f'{x:.4f}%'),
                textposition='outside',
                hovertemplate=(
                    '<b>%{y}</b><br>' +
                    'Abundance: %{x:.4f}%<br>' +
                    'Risk: %{customdata[0]}<br>' +
                    'Confidence: %{customdata[1]:.3f}<extra></extra>'
                ),
                customdata=list(zip(df['risk_level'], df['confidence']))
            ))
            
            fig.update_layout(
                title="Top Pathogenic Organisms by Relative Abundance<br><sub>Color-coded by risk level</sub>",
                xaxis_title="Relative Abundance (%)",
                yaxis_title="Organism",
                template="plotly_white",
                height=max(500, len(df) * 35),
                margin=dict(l=280, r=120, t=100, b=60),
                font=dict(size=10),
                showlegend=False
            )
            
            filepath = self.save_plot(fig, "pathogen_abundance_heatmap.html")
            self.formatter.success(f"Abundance heatmap saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.warning(f"Error creating abundance heatmap: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_all_visualizations(self, report_path):
        """Create all pathogen visualizations"""
        self.formatter.section_header("Creating Pathogen Visualizations")
        
        viz_files = []
        
        risk_file = self.create_risk_assessment_chart(report_path)
        if risk_file:
            viz_files.append(risk_file)
        
        who_file = self.create_who_priority_distribution(report_path)
        if who_file:
            viz_files.append(who_file)
        
        conf_file = self.create_detection_confidence_chart(report_path)
        if conf_file:
            viz_files.append(conf_file)
        
        abund_file = self.create_abundance_heatmap(report_path)
        if abund_file:
            viz_files.append(abund_file)
        
        self.formatter.success(f"Generated {len(viz_files)} pathogen visualizations")
        return viz_files