from pathlib import Path
from datetime import datetime
from .base_visualizer import BaseVisualizer

class DashboardGenerator(BaseVisualizer):
    """Generates the final HTML analysis dashboard with a sleek and professional design."""
    
    def __init__(self, output_dir: Path, analysis_type: str):
        super().__init__(output_dir)
        self.analysis_type = analysis_type
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.files_exist = self._check_files()
        self.theme = {
            'color': "#4CAF50" if self.analysis_type == 'fastq' else "#007bff",
            'name': "FASTQ Analysis" if self.analysis_type == 'fastq' else "FASTA Analysis",
            'version': "v3.6.0"
        }

    def _check_files(self) -> dict:
        """Checks for the existence of all key report files to display on the dashboard."""
        files_dict = {
            # Common Taxonomic Reports
            'taxonomic_report': (self.output_dir / "taxonomic_classification_report.txt").exists(),
            'taxonomic_abundance_chart': (self.output_dir / "taxonomic_abundance_chart.html").exists(),
            'taxonomy_krona': (self.output_dir / "taxonomy_krona.html").exists(),
            
            # Common Functional Reports
            'functional_report': (self.output_dir / "functional_annotation_report.txt").exists(),
            'functional_annotations_tsv': (self.output_dir / "functional_annotations.tsv").exists(),
            'protein_length_analysis': (self.output_dir / "protein_length_analysis.html").exists(),
            'functional_categories': (self.output_dir / "functional_categories.html").exists(),
            
            # Common Pathogen Detection Reports
            'pathogen_detection_report': (self.output_dir / "pathogen_detection_report.txt").exists(),
            
            # New Pathogen Visualizations (Plotly HTML)
            'pathogen_risk_assessment': (self.output_dir / "pathogen_risk_assessment.html").exists(),
            'who_priority_distribution': (self.output_dir / "who_priority_distribution.html").exists(),
            'detection_confidence': (self.output_dir / "detection_confidence.html").exists(),
            'diversity_metrics': (self.output_dir / "diversity_metrics.html").exists(),
        }
        
        if self.analysis_type == 'fastq':
            # FASTQ-specific files (Kraken2/Bracken)
            files_dict.update({
                'bracken_report_txt': (self.output_dir / "bracken_report.txt").exists(),
                'bracken_report_tsv': (self.output_dir / "bracken_report.tsv").exists(),
                'kraken_report': (self.output_dir / "kraken_report.txt").exists(),
                'kraken_classified': (self.output_dir / "kraken_classified.txt").exists(),
                'pathogen_results': (self.output_dir / "pathogen_results.txt").exists(),
            })
        elif self.analysis_type == 'fasta':
            # FASTA-specific files (BLAST + ML)
            files_dict.update({
                'blast_ml_pathogen_report': (self.output_dir / "blast_ml_pathogen_report.txt").exists(),
                'ml_predictions_csv': (self.output_dir / "ml_pathogen_predictions.csv").exists(),
            })
        
        return files_dict

    def _get_styles(self) -> str:
        """Returns the CSS styles for the dashboard."""
        return f"""
        <style>
            :root {{
                --primary-color: {self.theme['color']};
                --primary-light: {self.theme['color']}15;
                --primary-dark: {self.theme['color']}dd;
                --text-primary: #1a202c;
                --text-secondary: #4a5568;
                --text-muted: #718096;
                --background: #f7fafc;
                --card-background: #ffffff;
                --border: #e2e8f0;
                --border-light: #f1f5f9;
                --success: #38a169;
                --warning: #d69e2e;
                --danger: #e53e3e;
                --info: #3182ce;
                --shadow-sm: 0 1px 3px rgba(0, 0, 0, 0.1), 0 1px 2px rgba(0, 0, 0, 0.06);
                --shadow-md: 0 4px 6px rgba(0, 0, 0, 0.07), 0 2px 4px rgba(0, 0, 0, 0.06);
                --shadow-lg: 0 10px 15px rgba(0, 0, 0, 0.1), 0 4px 6px rgba(0, 0, 0, 0.05);
                --shadow-xl: 0 20px 25px rgba(0, 0, 0, 0.1), 0 10px 10px rgba(0, 0, 0, 0.04);
            }}

            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}

            body {{
                font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
                line-height: 1.6;
                color: var(--text-primary);
                background: #ffffff;
                min-height: 100vh;
                -webkit-font-smoothing: antialiased;
                -moz-osx-font-smoothing: grayscale;
            }}

            .container {{
                max-width: 1400px;
                margin: 0 auto;
                padding: 2rem 1rem;
            }}

            .header {{
                background: rgba(255, 255, 255, 0.95);
                backdrop-filter: blur(20px);
                border-radius: 20px;
                padding: 3rem 2rem;
                margin-bottom: 2.5rem;
                box-shadow: var(--shadow-xl);
                text-align: center;
                position: relative;
                overflow: hidden;
                border: 1px solid rgba(255, 255, 255, 0.2);
            }}

            .header::before {{
                content: '';
                position: absolute;
                top: 0;
                left: 0;
                right: 0;
                height: 5px;
                background: linear-gradient(90deg, var(--primary-color), var(--primary-dark));
            }}

            .header::after {{
                content: '';
                position: absolute;
                top: -50%;
                right: -50%;
                width: 200%;
                height: 200%;
                background: radial-gradient(circle, var(--primary-light) 0%, transparent 50%);
                opacity: 0.3;
                pointer-events: none;
            }}

            .header-content {{
                position: relative;
                z-index: 2;
            }}

            .header h1 {{
                font-size: 3rem;
                font-weight: 800;
                margin-bottom: 0.5rem;
                background: linear-gradient(135deg, var(--primary-color), var(--primary-dark));
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                background-clip: text;
                letter-spacing: -0.02em;
            }}

            .analysis-type {{
                display: inline-block;
                background: var(--primary-color);
                color: white;
                padding: 0.75rem 1.5rem;
                border-radius: 50px;
                font-weight: 600;
                font-size: 0.95rem;
                margin-bottom: 1.5rem;
                box-shadow: var(--shadow-md);
                text-transform: uppercase;
                letter-spacing: 0.05em;
            }}

            .header-meta {{
                display: flex;
                justify-content: center;
                gap: 2rem;
                margin-top: 1rem;
                flex-wrap: wrap;
            }}

            .header-meta span {{
                background: var(--border-light);
                color: var(--text-secondary);
                padding: 0.5rem 1rem;
                border-radius: 8px;
                font-size: 0.9rem;
                font-weight: 500;
                display: flex;
                align-items: center;
                gap: 0.5rem;
            }}

            .section {{
                margin-bottom: 2.5rem;
            }}

            .section-header {{
                background: rgba(255, 255, 255, 0.9);
                backdrop-filter: blur(10px);
                border-radius: 16px;
                padding: 1.25rem 1.5rem;
                margin-bottom: 1.5rem;
                border-left: 5px solid var(--primary-color);
                box-shadow: var(--shadow-md);
                position: relative;
            }}

            .section h2 {{
                font-size: 1.5rem;
                font-weight: 700;
                color: var(--text-primary);
                margin: 0;
                display: flex;
                align-items: center;
                gap: 0.75rem;
            }}

            .section-subtitle {{
                color: var(--text-secondary);
                font-size: 0.9rem;
                margin-top: 0.4rem;
                font-weight: 500;
            }}

            .grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                gap: 1.25rem;
                margin-bottom: 1rem;
            }}

            .card {{
                background: rgba(255, 255, 255, 0.95);
                backdrop-filter: blur(10px);
                border-radius: 12px;
                padding: 1.25rem;
                box-shadow: var(--shadow-md);
                transition: all 0.3s cubic-bezier(0.175, 0.885, 0.32, 1.275);
                border: 1px solid rgba(255, 255, 255, 0.2);
                position: relative;
                overflow: hidden;
                min-height: 140px;
                display: flex;
                flex-direction: column;
                justify-content: space-between;
            }}

            .card::before {{
                content: '';
                position: absolute;
                top: 0;
                left: 0;
                right: 0;
                height: 3px;
                background: linear-gradient(90deg, var(--primary-color), var(--primary-dark));
                transform: scaleX(0);
                transition: transform 0.4s ease;
                transform-origin: left;
            }}

            .card::after {{
                content: '';
                position: absolute;
                top: 0;
                right: 0;
                width: 80px;
                height: 80px;
                background: radial-gradient(circle, var(--primary-light) 0%, transparent 70%);
                transform: translate(25px, -25px);
                opacity: 0;
                transition: opacity 0.3s ease;
            }}

            .card:hover {{
                transform: translateY(-4px) scale(1.01);
                box-shadow: var(--shadow-lg);
            }}

            .card:hover::before {{
                transform: scaleX(1);
            }}

            .card:hover::after {{
                opacity: 1;
            }}

            .card-content {{
                position: relative;
                z-index: 2;
                flex-grow: 1;
            }}

            .card h3 {{
                font-size: 1.1rem;
                font-weight: 700;
                color: var(--text-primary);
                margin-bottom: 0.75rem;
                display: flex;
                align-items: center;
                gap: 0.6rem;
                line-height: 1.3;
            }}

            .icon {{
                font-size: 1.3rem;
                background: var(--primary-light);
                padding: 0.4rem;
                border-radius: 8px;
                display: flex;
                align-items: center;
                justify-content: center;
                min-width: 2.4rem;
                height: 2.4rem;
            }}

            .card p {{
                color: var(--text-secondary);
                margin-bottom: 1rem;
                font-size: 0.85rem;
                line-height: 1.5;
                flex-grow: 1;
            }}

            .card-footer {{
                display: flex;
                justify-content: space-between;
                align-items: center;
                margin-top: auto;
                padding-top: 0.75rem;
                border-top: 1px solid var(--border-light);
            }}

            .card a {{
                display: inline-flex;
                align-items: center;
                gap: 0.4rem;
                color: white;
                background: var(--primary-color);
                text-decoration: none;
                font-weight: 600;
                padding: 0.6rem 1.2rem;
                border-radius: 8px;
                transition: all 0.3s ease;
                font-size: 0.8rem;
                box-shadow: var(--shadow-sm);
                text-transform: uppercase;
                letter-spacing: 0.025em;
            }}

            .card a:hover {{
                background: var(--primary-dark);
                transform: translateX(2px);
                box-shadow: var(--shadow-md);
            }}

            .file-type {{
                background: var(--border-light);
                color: var(--text-muted);
                padding: 0.2rem 0.6rem;
                border-radius: 16px;
                font-size: 0.75rem;
                font-weight: 500;
                text-transform: uppercase;
            }}

            .file-type.interactive {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
            }}

            .file-type.static {{
                background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
                color: white;
            }}

            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 1rem;
                margin: 1.5rem 0;
                padding: 1.25rem;
                background: rgba(255, 255, 255, 0.9);
                border-radius: 16px;
                backdrop-filter: blur(10px);
            }}

            .stat-item {{
                text-align: center;
                padding: 0.75rem;
                border-radius: 12px;
                background: var(--border-light);
            }}

            .stat-number {{
                font-size: 1.3rem;
                font-weight: 700;
                color: var(--primary-color);
                display: block;
            }}

            .stat-label {{
                font-size: 0.8rem;
                color: var(--text-secondary);
                margin-top: 0.25rem;
                text-transform: uppercase;
                letter-spacing: 0.05em;
            }}

            .badge {{
                display: inline-block;
                padding: 0.25rem 0.6rem;
                border-radius: 12px;
                font-size: 0.7rem;
                font-weight: 600;
                margin-left: 0.5rem;
                text-transform: uppercase;
                letter-spacing: 0.03em;
            }}

            .badge.new {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
            }}

            .badge.enhanced {{
                background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
                color: white;
            }}

            .footer {{
                background: rgba(255, 255, 255, 0.9);
                backdrop-filter: blur(10px);
                text-align: center;
                padding: 1.5rem;
                border-radius: 16px;
                color: var(--text-secondary);
                font-size: 0.9rem;
                margin-top: 2.5rem;
                box-shadow: var(--shadow-md);
            }}

            .footer-logo {{
                font-weight: 700;
                color: var(--primary-color);
                margin-bottom: 0.5rem;
            }}

            @media (max-width: 768px) {{
                .container {{
                    padding: 1rem;
                }}

                .header {{
                    padding: 2rem 1.5rem;
                }}

                .header h1 {{
                    font-size: 2.5rem;
                }}

                .header-meta {{
                    flex-direction: column;
                    gap: 1rem;
                }}

                .grid {{
                    grid-template-columns: 1fr;
                }}

                .card {{
                    padding: 1rem;
                    min-height: 120px;
                }}

                .section-header {{
                    padding: 1rem;
                }}

                .section h2 {{
                    font-size: 1.3rem;
                }}
            }}

            /* Animations */
            @keyframes slideInUp {{
                from {{
                    opacity: 0;
                    transform: translateY(20px);
                }}
                to {{
                    opacity: 1;
                    transform: translateY(0);
                }}
            }}

            .card {{
                animation: slideInUp 0.5s ease forwards;
            }}

            .card:nth-child(2) {{ animation-delay: 0.1s; }}
            .card:nth-child(3) {{ animation-delay: 0.15s; }}
            .card:nth-child(4) {{ animation-delay: 0.2s; }}
            .card:nth-child(5) {{ animation-delay: 0.25s; }}
            .card:nth-child(6) {{ animation-delay: 0.3s; }}
            .card:nth-child(7) {{ animation-delay: 0.35s; }}
            .card:nth-child(8) {{ animation-delay: 0.4s; }}
        </style>
        """

    def _create_card(self, icon: str, title: str, description: str, link_key: str, 
                    file_type: str = "TXT", badge: str = None) -> str:
        """Creates the HTML for a single card if the corresponding file exists."""
        if self.files_exist.get(link_key, False):
            filename_map = {
                # Taxonomic Reports
                'taxonomic_report': 'taxonomic_classification_report.txt',
                'taxonomic_abundance_chart': 'taxonomic_abundance_chart.html',
                'taxonomy_krona': 'taxonomy_krona.html',
                
                # Functional Reports
                'functional_report': 'functional_annotation_report.txt',
                'functional_annotations_tsv': 'functional_annotations.tsv',
                'protein_length_analysis': 'protein_length_analysis.html',
                'functional_categories': 'functional_categories.html',
                
                # Pathogen Detection Reports
                'pathogen_detection_report': 'pathogen_detection_report.txt',
                
                # New Interactive Pathogen Visualizations (Plotly)
                'pathogen_risk_assessment': 'pathogen_risk_assessment.html',
                'who_priority_distribution': 'who_priority_distribution.html',
                'detection_confidence': 'detection_confidence.html',
                'diversity_metrics': 'diversity_metrics.html',
                
                # FASTQ-specific
                'bracken_report_txt': 'bracken_report.txt',
                'bracken_report_tsv': 'bracken_report.tsv',
                'kraken_report': 'kraken_report.txt',
                'kraken_classified': 'kraken_classified.txt',
                'pathogen_results': 'pathogen_results.txt',
                
                # FASTA-specific
                'blast_ml_pathogen_report': 'blast_ml_pathogen_report.txt',
                'ml_predictions_csv': 'ml_pathogen_predictions.csv',
            }
            
            filename = filename_map.get(link_key, '')
            
            # Determine file type badge style
            file_type_class = ""
            if file_type == "HTML":
                file_type_class = "interactive"
            elif file_type == "PNG":
                file_type_class = "static"
            
            badge_html = ""
            if badge:
                badge_class = badge.lower()
                badge_html = f'<span class="badge {badge_class}">{badge}</span>'

            return f"""
            <div class="card">
                <div class="card-content">
                    <h3><span class="icon">{icon}</span> {title}{badge_html}</h3>
                    <p>{description}</p>
                </div>
                <div class="card-footer">
                    <span class="file-type {file_type_class}">{file_type}</span>
                    <a href="{filename}" target="_blank">Open Report &rarr;</a>
                </div>
            </div>
            """
        return ""

    def create_dashboard(self):
        """Generates and saves the complete HTML dashboard."""
        html_head = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>MetaQuest Analysis Dashboard - {self.theme['name']}</title>
            <link rel="preconnect" href="https://fonts.googleapis.com">
            <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
            <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet">
            {self._get_styles()}
        </head>
        <body>
            <div class="container">
        """
        
        html_header_body = self._create_header_and_body()
        
        html_foot = f"""
                <div class="footer">
                    <div class="footer-logo">🧬 MetaQuest Bioinformatics Platform</div>
                    <p>Advanced metagenomics analysis pipeline combining traditional microbiology with cutting-edge AI technologies</p>
                    <p>Version {self.theme['version']} • Powered by Kraken2, BLAST, Prokka & Custom ML Models</p>
                    <p style="margin-top: 0.5rem; font-size: 0.85rem;">
                        Interactive visualizations powered by Plotly | Static charts generated with Matplotlib
                    </p>
                </div>
            </div>
        </body>
        </html>
        """
        
        final_html = html_head + html_header_body + html_foot
        
        dashboard_path = self.output_dir / "analysis_dashboard.html"
        with open(dashboard_path, 'w', encoding='utf-8') as f:
            f.write(final_html)
        print(f"✓ Analysis dashboard created at: {dashboard_path}")

    def _create_header_and_body(self) -> str:
        """Creates the header and main body content based on analysis type."""
        # Count available files for stats
        total_files = sum(1 for exists in self.files_exist.values() if exists)
        report_files = sum(1 for k, v in self.files_exist.items() if v and ('report' in k or 'summary' in k))
        viz_files = sum(1 for k, v in self.files_exist.items() if v and ('chart' in k or 'distribution' in k or 'krona' in k))
        
        header = f"""
        <div class="header">
            <div class="header-content">
                <h1>🧬 MetaQuest Analysis Report</h1>
                <div class="analysis-type">{self.theme['name']}</div>
                <div class="header-meta">
                    <span>📅 {self.timestamp}</span>
                    <span>🔬 Version {self.theme['version']}</span>
                    <span>📊 {total_files} Files Generated</span>
                </div>
            </div>
        </div>
        
        <div class="stats-grid">
            <div class="stat-item">
                <span class="stat-number">{report_files}</span>
                <span class="stat-label">Reports</span>
            </div>
            <div class="stat-item">
                <span class="stat-number">{viz_files}</span>
                <span class="stat-label">Visualizations</span>
            </div>
            <div class="stat-item">
                <span class="stat-number">{self.theme['version']}</span>
                <span class="stat-label">Pipeline Version</span>
            </div>
        </div>
        """
        
        taxonomic_section, functional_section, pathogen_section = "", "", ""

        if self.analysis_type == 'fastq':
            taxonomic_section = self._create_fastq_taxonomic_section()
            functional_section = self._create_functional_section()
            pathogen_section = self._create_fastq_pathogen_section()
            
        elif self.analysis_type == 'fasta':
            taxonomic_section = self._create_fasta_taxonomic_section()
            functional_section = self._create_functional_section()
            pathogen_section = self._create_fasta_pathogen_section()
        
        return header + taxonomic_section + functional_section + pathogen_section

    def _create_fastq_taxonomic_section(self) -> str:
        """Create taxonomic section for FASTQ analysis"""
        return f"""
        <div class="section">
            <div class="section-header">
                <h2>📊 Taxonomic Classification</h2>
                <p class="section-subtitle">Species identification using Kraken2/Bracken k-mer profiling</p>
            </div>
            <div class="grid">
                {self._create_card("📋", "Taxonomic Classification Report", "Comprehensive species abundance and classification results.", "taxonomic_report", "TXT")}
                {self._create_card("📈", "Abundance Chart", "Interactive visualization of taxonomic abundance distribution.", "taxonomic_abundance_chart", "HTML")}
                {self._create_card("🌍", "Krona Taxonomy Chart", "Hierarchical interactive taxonomy browser.", "taxonomy_krona", "HTML")}
                {self._create_card("📊", "Bracken Report (Text)", "Human-readable Bracken species abundance report.", "bracken_report_txt", "TXT")}
                {self._create_card("📄", "Bracken Report (Data)", "Machine-readable species-level abundance data.", "bracken_report_tsv", "TSV")}
                {self._create_card("📑", "Kraken Classification Report", "Raw Kraken2 classification results.", "kraken_report", "TXT")}
            </div>
        </div>
        """

    def _create_fasta_taxonomic_section(self) -> str:
        """Create taxonomic section for FASTA analysis"""
        return f"""
        <div class="section">
            <div class="section-header">
                <h2>📊 Taxonomic Classification</h2>
                <p class="section-subtitle">BLAST-based organism identification and abundance analysis</p>
            </div>
            <div class="grid">
                {self._create_card("📋", "Taxonomic Classification Report", "BLAST-based species identification results.", "taxonomic_report", "TXT")}
                {self._create_card("📈", "Abundance Chart", "Interactive organism abundance distribution.", "taxonomic_abundance_chart", "HTML")}
                {self._create_card("🌍", "Krona Taxonomy Chart", "Hierarchical taxonomy visualization.", "taxonomy_krona", "HTML")}
            </div>
        </div>
        """

    def _create_functional_section(self) -> str:
        """Create functional annotation section (common to both)"""
        return f"""
        <div class="section">
            <div class="section-header">
                <h2>🧪 Functional Annotation</h2>
                <p class="section-subtitle">Gene prediction and protein functional characterization</p>
            </div>
            <div class="grid">
                {self._create_card("🔬", "Functional Annotation Report", "Comprehensive gene and protein annotation summary.", "functional_report", "TXT")}
                {self._create_card("📊", "Functional Annotations (Data)", "Detailed annotation data for all predicted features.", "functional_annotations_tsv", "TSV")}
                {self._create_card("📏", "Protein Length Analysis", "Statistical analysis of predicted protein sizes.", "protein_length_analysis", "HTML")}
                {self._create_card("📈", "Functional Categories", "Distribution of protein functional classifications.", "functional_categories", "HTML")}
            </div>
        </div>
        """

    def _create_fastq_pathogen_section(self) -> str:
        """Create pathogen section for FASTQ analysis"""
        return f"""
        <div class="section">
            <div class="section-header">
                <h2>🦠 Pathogen Detection & Risk Assessment</h2>
                <p class="section-subtitle">Multi-method pathogen screening with clinical risk stratification</p>
            </div>
            <div class="grid">
                {self._create_card("📋", "Pathogen Detection Report", "Comprehensive pathogen analysis with WHO priority classification.", "pathogen_detection_report", "TXT")}
                {self._create_card("🧬", "Pathogen Results", "Detailed pathogen detection results from multiple methods.", "pathogen_results", "TXT")}
                {self._create_card("📑", "Kraken Classified Sequences", "Raw classified sequence data from Kraken2.", "kraken_classified", "TXT")}
            </div>
            
            <div class="section-header" style="margin-top: 2rem;">
                <h2>📊 Interactive Visualizations <span class="badge new">Enhanced</span></h2>
                <p class="section-subtitle">Interactive Plotly charts with advanced filtering and zoom capabilities</p>
            </div>
            <div class="grid">
                {self._create_card("🎯", "Risk Assessment Chart", "Interactive pathogen risk level analysis and classification.", "pathogen_risk_assessment", "HTML", "NEW")}
                {self._create_card("🌍", "WHO Priority Distribution", "WHO 2024 priority pathogen classification breakdown.", "who_priority_distribution", "HTML", "NEW")}
                {self._create_card("📈", "Detection Confidence", "Confidence scores and multi-method validation analysis.", "detection_confidence", "HTML", "NEW")}
                {self._create_card("📊", "Diversity Metrics", "Microbial diversity and richness analysis.", "diversity_metrics", "HTML", "NEW")}
            </div>
        </div>
        """

    def _create_fasta_pathogen_section(self) -> str:
        """Create pathogen section for FASTA analysis"""
        return f"""
        <div class="section">
            <div class="section-header">
                <h2>🤖 Integrated Pathogen Analysis (BLAST + ML)</h2>
                <p class="section-subtitle">Dual-method pathogen detection combining sequence alignment and machine learning</p>
            </div>
            <div class="grid">
                {self._create_card("📋", "BLAST+ML Pathogen Report", "Integrated pathogen analysis combining BLAST taxonomy and ML predictions.", "blast_ml_pathogen_report", "TXT")}
                {self._create_card("🤖", "ML Predictions (CSV)", "Individual protein pathogenicity predictions with confidence scores.", "ml_predictions_csv", "CSV")}
            </div>
            
            <div class="section-header" style="margin-top: 2rem;">
                <h2>📊 Interactive Visualizations <span class="badge new">Enhanced</span></h2>
                <p class="section-subtitle">Interactive Plotly charts with drill-down capabilities</p>
            </div>
            <div class="grid">
                {self._create_card("🎯", "Risk Assessment Chart", "Integrated BLAST+ML pathogen risk classification.", "pathogen_risk_assessment", "HTML", "NEW")}
                {self._create_card("🌍", "WHO Priority Distribution", "WHO priority pathogen breakdown from BLAST results.", "who_priority_distribution", "HTML", "NEW")}
                {self._create_card("📈", "Detection Confidence", "ML prediction confidence and BLAST quality metrics.", "detection_confidence", "HTML", "NEW")}
                {self._create_card("📊", "Diversity Analysis", "Taxonomic diversity and community structure metrics.", "diversity_metrics", "HTML", "NEW")}
            </div>
        </div>
                <h2>📸 Static Charts</h2>
                <p class="section-subtitle">Publication-ready static visualizations</p>
            </div>
            <div class="grid">
                {self._create_card("⚠️", "Risk Distribution (PNG)", "Static pathogen risk level distribution.", "pathogen_risk_distribution", "PNG")}
                {self._create_card("📊", "Top Pathogens (PNG)", "Static bar chart of detected pathogenic organisms.", "pathogen_abundance_top15", "PNG")}
                {self._create_card("🏥", "WHO Classification (PNG)", "Static WHO priority classification chart.", "pathogen_who_classification", "PNG")}
                {self._create_card("🎲", "ML Confidence (PNG)", "Static ML prediction confidence distribution.", "pathogen_confidence_methods", "PNG")}
            </div>
        </div>
        """


# --- Main wrapper function ---
def create_dashboard(analysis_type: str, output_dir: Path):
    """Main entry point to create the final analysis dashboard."""
    generator = DashboardGenerator(output_dir, analysis_type)
    generator.create_dashboard()
    return str(output_dir / "analysis_dashboard.html")