from pathlib import Path
from datetime import datetime
from .visualization import BaseVisualizer

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
            'version': "v3.2.2" # You can make this dynamic if needed
        }

    def _check_files(self) -> dict:
        """Checks for the existence of all key report files to display on the dashboard."""
        files_dict = {
            # Common files for both FASTA and FASTQ
            'taxonomic_report': (self.output_dir / "taxonomic_classification_report.txt").exists(),
            'functional_report': (self.output_dir / "functional_annotation_report.txt").exists(),
            'protein_length_distribution': (self.output_dir / "protein_length_distribution.html").exists(),
            'pathogen_risk_chart': (self.output_dir / "pathogen_risk_detection.html").exists(),
            'taxonomic_abundance_chart': (self.output_dir / "taxonomic_abundance_chart.html").exists(),
            'taxonomy_krona': (self.output_dir / "taxonomy_krona.html").exists(),
            'detection_method_coverage': (self.output_dir / "detection_method_coverage.html").exists(),
        }
        
        if self.analysis_type == 'fastq':
            # FASTQ-specific files (traditional pathogen screening)
            files_dict.update({
                'pathogen_summary': (self.output_dir / "pathogen_summary.txt").exists(),
                'bracken_report': (self.output_dir / "bracken_report.tsv").exists(),
                'kraken_report': (self.output_dir / "kraken_report.txt").exists(),
                'amr_hits': (self.output_dir / "amr_hits.tsv").exists(),
                'pathogen_results': (self.output_dir / "pathogen_results.txt").exists(),
                'Pathogenic_Proteins': (self.output_dir / "Pathogenic_Proteins.txt").exists(),
            })
        elif self.analysis_type == 'fasta':
            # FASTA-specific files (BLAST + ML)
            files_dict.update({
                'blast_ml_summary': (self.output_dir / "blast_ml_pathogen_summary.txt").exists(),
                'ml_predictions_csv': (self.output_dir / "ml_pathogen_predictions.csv").exists(),
                'organism_comparison_csv': (self.output_dir / "organism_comparison_data.csv").exists(),
                'blast_ml_integrated_report': (self.output_dir / "blast_ml_integrated_pathogen_report.json").exists(),
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
        </style>
        """

    def _create_card(self, icon: str, title: str, description: str, link_key: str, file_type: str = "TXT") -> str:
        """Creates the HTML for a single card if the corresponding file exists."""
        if self.files_exist.get(link_key, False):
            filename_map = {
                # Common files
                'taxonomic_report': 'taxonomic_classification_report.txt',
                'functional_report': 'functional_annotation_report.txt',
                'protein_length_distribution': 'protein_length_distribution.html',
                'pathogen_risk_chart': 'pathogen_risk_detection.html',
                'taxonomic_abundance_chart': 'taxonomic_abundance_chart.html',
                'taxonomy_krona': 'taxonomy_krona.html',
                'detection_method_coverage': 'detection_method_coverage.html',
                
                # FASTQ-specific files
                'pathogen_summary': 'pathogen_summary.txt',
                'bracken_report': 'bracken_report.tsv',
                'pathogen_results': 'pathogen_results.txt',
                'Pathogenic_Proteins': 'Pathogenic_Proteins.txt',
                
                # FASTA-specific files
                'blast_ml_summary': 'blast_ml_pathogen_summary.txt',
                'ml_predictions_csv': 'ml_pathogen_predictions.csv',
                'organism_comparison_csv': 'organism_comparison_data.csv',
                'blast_ml_integrated_report': 'blast_ml_integrated_pathogen_report.json',
            }
            
            filename = filename_map.get(link_key, '')

            return f"""
            <div class="card">
                <div class="card-content">
                    <h3><span class="icon">{icon}</span> {title}</h3>
                    <p>{description}</p>
                </div>
                <div class="card-footer">
                    <span class="file-type">{file_type}</span>
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
            <title>MetaQuest Analysis Dashboard</title>
            <link rel="preconnect" href="https://fonts.googleapis.com">
            <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
            <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
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
                </div>
            </div>
        </body>
        </html>
        """
        
        final_html = html_head + html_header_body + html_foot
        
        with open(self.output_dir / "analysis_dashboard.html", 'w') as f:
            f.write(final_html)
        print(f"✓ Analysis dashboard created at: {self.output_dir / 'analysis_dashboard.html'}")

    def _create_header_and_body(self) -> str:
        """Creates the header and main body content based on analysis type."""
        header = f"""
        <div class="header">
            <h1>🧬 MetaQuest Report</h1>
            <div class="analysis-type">{self.theme['name']}</div>
            <p>Version: {self.theme['version']} | Generated on: {self.timestamp}</p>
        </div>
        """
        
        taxonomic_section, functional_section, pathogen_section = "", "", ""

        if self.analysis_type == 'fastq':
            taxonomic_section = f"""
            <div class="section">
                <div class="section-header">
                    <h2>📊 Taxonomic Analysis</h2>
                </div>
                <div class="grid">
                    {self._create_card("📋", "Taxonomic Classification Report", "Comprehensive species identification and abundance analysis.", "taxonomic_report", "TXT")}
                    {self._create_card("📈", "Taxonomic Abundance Chart", "Interactive visualization of the most abundant taxa in your sample.", "taxonomic_abundance_chart", "HTML")}
                    {self._create_card("🌐", "Krona Taxonomy Chart", "Hierarchical interactive view of taxonomic classification.", "taxonomy_krona", "HTML")}
                    {self._create_card("📊", "Bracken Report", "Refined abundance estimates at species level from Bracken.", "bracken_report", "TSV")}
                </div>
            </div>
            """
            
            functional_section = f"""
            <div class="section">
                <div class="section-header">
                    <h2>🧪 Functional Analysis</h2>
                </div>
                <div class="grid">
                    {self._create_card("🔬", "Functional Annotation Report", "Gene prediction and protein annotation summary from Prokka.", "functional_report", "TXT")}
                    {self._create_card("📏", "Protein Length Distribution", "Statistical distribution of predicted protein lengths.", "protein_length_distribution", "HTML")}
                </div>
            </div>
            """
            
            pathogen_section = f"""
            <div class="section">
                <div class="section-header">
                    <h2>🦠 Pathogen Analysis (Traditional Screening)</h2>
                </div>
                <div class="grid">
                    {self._create_card("🎯", "Pathogen Summary Report", "Main clinical pathogen detection summary with risk assessment.", "pathogen_summary", "TXT")}
                    {self._create_card("☣️", "Pathogen Risk Detection Chart", "Interactive visualization of detected pathogens by risk level.", "pathogen_risk_chart", "HTML")}
                    {self._create_card("📄", "Pathogenic Proteins", "Comprehensive pathogen detection results and annotations.", "Pathogenic_Proteins", "TXT")}
                    {self._create_card("📈", "Detection Method Coverage", "Coverage analysis of different pathogen detection methods.", "detection_method_coverage", "HTML")}
                </div>
            </div>
            """
            
        elif self.analysis_type == 'fasta':
            taxonomic_section = f"""
            <div class="section">
                <div class="section-header">
                    <h2>📊 Taxonomic Analysis (BLAST-based)</h2>
                </div>
                <div class="grid">
                    {self._create_card("📋", "Taxonomic Classification Report", "Species identification results from BLAST analysis.", "taxonomic_report", "TXT")}
                    {self._create_card("📈", "Taxonomic Abundance Chart", "Interactive chart showing organism abundance distribution.", "taxonomic_abundance_chart", "HTML")}
                    {self._create_card("🌐", "Krona Taxonomy Chart", "Hierarchical view of taxonomic classification results.", "taxonomy_krona", "HTML")}
                    {self._create_card("📊", "Organism Comparison Data", "Detailed BLAST hit statistics and organism comparison metrics.", "organism_comparison_csv", "CSV")}
                </div>
            </div>
            """
            
            functional_section = f"""
            <div class="section">
                <div class="section-header">
                    <h2>🧪 Functional Analysis</h2>
                </div>
                <div class="grid">
                    {self._create_card("🔬", "Functional Annotation Report", "Gene prediction and protein annotation summary from Prokka.", "functional_report", "TXT")}
                    {self._create_card("📏", "Protein Length Distribution", "Statistical distribution of predicted protein lengths.", "protein_length_distribution", "HTML")}
                </div>
            </div>
            """
            
            pathogen_section = f"""
            <div class="section">
                <div class="section-header">
                    <h2>🤖 Pathogen Analysis (BLAST + Machine Learning)</h2>
                </div>
                <div class="grid">
                    {self._create_card("🎯", "BLAST+ML Integrated Summary", "Comprehensive pathogen analysis combining BLAST results with ML predictions.", "blast_ml_summary", "TXT")}
                    {self._create_card("💻", "ML Pathogenicity Predictions", "Individual protein pathogenicity predictions from machine learning models.", "ml_predictions_csv", "CSV")}
                    {self._create_card("☣️", "Pathogen Risk Detection Chart", "Interactive visualization of detected pathogens and risk levels.", "pathogen_risk_chart", "HTML")}
                    {self._create_card("📈", "Detection Method Coverage", "Analysis of detection method coverage and effectiveness.", "detection_method_coverage", "HTML")}
                </div>
            </div>
            """
        
        return header + taxonomic_section + functional_section + pathogen_section

# --- Main wrapper function ---
def create_dashboard(analysis_type: str, output_dir: Path):
    """Main entry point to create the final analysis dashboard."""
    generator = DashboardGenerator(output_dir, analysis_type)
    generator.create_dashboard()