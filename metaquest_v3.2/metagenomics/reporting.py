import pandas as pd
import json
from datetime import datetime
from pathlib import Path
from .config import *
from .visualization import *

def create_analysis_dashboard(output_dir):
    """Create comprehensive analysis dashboard with tabs"""
    dashboard_file = output_dir/"analysis_dashboard.html"
    
    html_files = list(output_dir.glob("*.html"))
    html_files = [f for f in html_files if f.name != "analysis_dashboard.html"]
    
    # Group files by category
    taxonomy_files = [f for f in html_files if any(term in f.name for term in ['taxonomy', 'kraken', 'bracken', 'krona'])]
    functional_files = [f for f in html_files if any(term in f.name for term in ['swissprot', 'prokka', 'functional', '3d_annotation'])]
    pathogen_files = [f for f in html_files if any(term in f.name for term in ['pathogen', 'blast'])]
    quality_files = [f for f in html_files if any(term in f.name for term in ['length', 'gc', 'quality', 'summary'])]
    other_files = [f for f in html_files if f not in taxonomy_files + functional_files + pathogen_files + quality_files]
    
    dashboard_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Metagenomics Analysis Dashboard</title>
        <style>
            body {{ 
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
                margin: 0; 
                padding: 20px; 
                background: #f8f9fa;
            }}
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 30px;
                border-radius: 10px;
                margin-bottom: 30px;
                text-align: center;
            }}
            .tab-container {{ 
                display: flex; 
                background: #fff; 
                border-radius: 10px 10px 0 0;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            .tab {{ 
                background: inherit; 
                border: none; 
                padding: 15px 25px; 
                cursor: pointer; 
                transition: all 0.3s;
                font-weight: 500;
            }}
            .tab:hover {{ background: #e9ecef; }}
            .tab.active {{ 
                background: #007bff; 
                color: white; 
                border-radius: 10px 10px 0 0;
            }}
            .tabcontent {{ 
                display: none; 
                padding: 30px; 
                background: white;
                border-radius: 0 0 10px 10px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                min-height: 400px;
            }}
            .active-tab {{ display: block; }}
            .grid {{ 
                display: grid; 
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
                gap: 20px; 
                margin-top: 20px;
            }}
            .card {{ 
                background: #f8f9fa; 
                padding: 20px; 
                border-radius: 8px; 
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                transition: transform 0.2s;
            }}
            .card:hover {{ transform: translateY(-2px); }}
            .card h3 {{ color: #495057; margin-top: 0; }}
            .card a {{ 
                color: #007bff; 
                text-decoration: none; 
                font-weight: bold;
                display: inline-block;
                margin-top: 10px;
                padding: 8px 16px;
                background: #e3f2fd;
                border-radius: 5px;
                transition: all 0.2s;
            }}
            .card a:hover {{ 
                background: #007bff; 
                color: white; 
            }}
            .summary-stats {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-bottom: 30px;
            }}
            .stat-card {{
                background: linear-gradient(135deg, #ff6b6b, #ee5a24);
                color: white;
                padding: 20px;
                border-radius: 10px;
                text-align: center;
            }}
            .stat-number {{ font-size: 2em; font-weight: bold; }}
            .stat-label {{ font-size: 0.9em; opacity: 0.9; }}
        </style>
        <script>
            function openTab(tabName) {{
                var i, tabcontent, tabs;
                tabcontent = document.getElementsByClassName("tabcontent");
                for (i = 0; i < tabcontent.length; i++) {{
                    tabcontent[i].classList.remove("active-tab");
                }}
                tabs = document.getElementsByClassName("tab");
                for (i = 0; i < tabs.length; i++) {{
                    tabs[i].classList.remove("active");
                }}
                document.getElementById(tabName).classList.add("active-tab");
                event.currentTarget.classList.add("active");
            }}
            
            // Auto-open first tab
            window.onload = function() {{
                document.querySelector('.tab').click();
            }}
        </script>
    </head>
    <body>
        <div class="header">
            <h1>🧬 Metagenomics Analysis Dashboard</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="tab-container">
            <button class="tab" onclick="openTab('summary')">📊 Summary</button>
            <button class="tab" onclick="openTab('taxonomy')">🦠 Taxonomy</button>
            <button class="tab" onclick="openTab('function')">⚙️ Function</button>
            <button class="tab" onclick="openTab('pathogen')">⚠️ Pathogens</button>
            <button class="tab" onclick="openTab('quality')">📈 Quality</button>
        </div>
        
        <div id="summary" class="tabcontent">
            <h2>Analysis Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <div class="stat-number">{len(html_files)}</div>
                    <div class="stat-label">Visualizations</div>
                </div>
                <div class="stat-card" style="background: linear-gradient(135deg, #4ecdc4, #26a69a);">
                    <div class="stat-number">{len(taxonomy_files)}</div>
                    <div class="stat-label">Taxonomy Plots</div>
                </div>
                <div class="stat-card" style="background: linear-gradient(135deg, #45b7d1, #2980b9);">
                    <div class="stat-number">{len(functional_files)}</div>
                    <div class="stat-label">Functional Plots</div>
                </div>
                <div class="stat-card" style="background: linear-gradient(135deg, #f39c12, #e67e22);">
                    <div class="stat-number">{len(pathogen_files)}</div>
                    <div class="stat-label">Pathogen Plots</div>
                </div>
            </div>
            
            <div class="grid">
                <div class="card">
                    <h3>🔍 Analysis Overview</h3>
                    <p>This dashboard provides comprehensive insights into your metagenomic sample, including taxonomic composition, functional annotations, and pathogen screening results.</p>
                </div>
                <div class="card">
                    <h3>📁 File Summary</h3>
                    <p>Total files generated: <strong>{len(html_files)}</strong></p>
                    <p>Navigate through the tabs above to explore different aspects of your analysis.</p>
                </div>
            </div>
        </div>
        
        <div id="taxonomy" class="tabcontent">
            <h2>🦠 Taxonomic Analysis</h2>
            <p>Explore the taxonomic composition of your sample</p>
            <div class="grid">
    """
    
    for html_file in taxonomy_files:
        dashboard_html += f"""
                <div class="card">
                    <h3>{html_file.stem.replace('_', ' ').title()}</h3>
                    <p>Interactive taxonomic visualization</p>
                    <a href="{html_file.name}" target="_blank">🔗 Open Visualization</a>
                </div>
        """
    
    dashboard_html += """
            </div>
        </div>
        
        <div id="function" class="tabcontent">
            <h2>⚙️ Functional Analysis</h2>
            <p>Functional annotation and protein analysis results</p>
            <div class="grid">
    """
    
    for html_file in functional_files:
        dashboard_html += f"""
                <div class="card">
                    <h3>{html_file.stem.replace('_', ' ').title()}</h3>
                    <p>Functional annotation visualization</p>
                    <a href="{html_file.name}" target="_blank">🔗 Open Visualization</a>
                </div>
        """
    
    dashboard_html += """
            </div>
        </div>
        
        <div id="pathogen" class="tabcontent">
            <h2>⚠️ Pathogen Screening</h2>
            <p>Pathogenic organism detection and analysis</p>
            <div class="grid">
    """
    
    for html_file in pathogen_files:
        dashboard_html += f"""
                <div class="card">
                    <h3>{html_file.stem.replace('_', ' ').title()}</h3>
                    <p>Pathogen screening results</p>
                    <a href="{html_file.name}" target="_blank">🔗 Open Visualization</a>
                </div>
        """
    
    dashboard_html += """
            </div>
        </div>
        
        <div id="quality" class="tabcontent">
            <h2>📈 Sequence Quality</h2>
            <p>Sequence statistics and quality metrics</p>
            <div class="grid">
    """
    
    for html_file in quality_files:
        dashboard_html += f"""
                <div class="card">
                    <h3>{html_file.stem.replace('_', ' ').title()}</h3>
                    <p>Quality and statistics visualization</p>
                    <a href="{html_file.name}" target="_blank">🔗 Open Visualization</a>
                </div>
        """
    
    dashboard_html += """
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(dashboard_file, 'w') as f:
        f.write(dashboard_html)
    print(f"✓ Created enhanced analysis dashboard: {dashboard_file}")

def generate_final_report(output_dir):
    """Generate comprehensive final report"""
    report_file = output_dir/"final_report.html"
    
    html_files = list(output_dir.glob("*.html"))
    txt_files = list(output_dir.glob("*.txt"))
    
    report_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Final Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .section {{ margin: 30px 0; padding: 20px; background: #f8f9fa; }}
        </style>
    </head>
    <body>
        <h1>Metagenomics Analysis Report</h1>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="section">
            <h2>Summary</h2>
            <p>Total visualizations: {len(html_files)}</p>
            <p>Data files: {len(txt_files)}</p>
        </div>
        
        <div class="section">
            <h2>Interactive Visualizations</h2>
            <ul>
    """
    
    for html_file in html_files:
        report_html += f'<li><a href="{html_file.name}" target="_blank">{html_file.stem}</a></li>\n'
    
    report_html += """
            </ul>
        </div>
    </body>
    </html>
    """
    
    with open(report_file, 'w') as f:
        f.write(report_html)
    print(f"✓ Generated final report: {report_file}") 