import subprocess
import os
import pandas as pd
from pathlib import Path
from .config import *
from .utils import check_dependencies

def generate_amr_report(amr_results, output_dir):
    """Generate antimicrobial resistance report from DIAMOND results"""
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    
    if not amr_results or not amr_results.exists():
        print("No AMR results to process")
        # Create empty report
        empty_report = {
            'summary': 'No antimicrobial resistance genes detected',
            'total_hits': 0,
            'resistance_classes': {},
            'top_hits': []
        }
        
        with open(output_dir / "antimicrobial_resistance_report.json", 'w') as f:
            import json
            json.dump(empty_report, f, indent=2)
        return
    
    try:
        # Read DIAMOND results
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle']
        df = pd.read_csv(amr_results, sep='\t', names=columns, header=None)
        
        if df.empty:
            print("No AMR hits found")
            return
        
        # Filter high-quality hits
        df_filtered = df[(df['pident'] >= 70) & (df['evalue'] <= 1e-10)]
        
        # Extract resistance information from titles
        resistance_classes = {}
        antibiotic_families = {
            'beta-lactam': ['beta-lactam', 'penicillin', 'ampicillin', 'cephalosporin'],
            'aminoglycoside': ['aminoglycoside', 'streptomycin', 'gentamicin', 'kanamycin'],
            'tetracycline': ['tetracycline', 'doxycycline'],
            'quinolone': ['quinolone', 'fluoroquinolone', 'ciprofloxacin'],
            'macrolide': ['macrolide', 'erythromycin', 'azithromycin'],
            'chloramphenicol': ['chloramphenicol'],
            'sulfonamide': ['sulfonamide', 'sulfamethoxazole'],
            'trimethoprim': ['trimethoprim'],
            'vancomycin': ['vancomycin'],
            'lincosamide': ['lincosamide', 'clindamycin']
        }
        
        for idx, row in df_filtered.iterrows():
            title_lower = row['stitle'].lower()
            for family, keywords in antibiotic_families.items():
                if any(keyword in title_lower for keyword in keywords):
                    if family not in resistance_classes:
                        resistance_classes[family] = []
                    resistance_classes[family].append({
                        'gene': row['sseqid'],
                        'identity': row['pident'],
                        'evalue': row['evalue'],
                        'description': row['stitle']
                    })
                    break
            else:
                # Unknown resistance mechanism
                if 'other' not in resistance_classes:
                    resistance_classes['other'] = []
                resistance_classes['other'].append({
                    'gene': row['sseqid'],
                    'identity': row['pident'],
                    'evalue': row['evalue'],
                    'description': row['stitle']
                })
        
        # Generate report
        amr_report = {
            'summary': f"Detected {len(df_filtered)} high-confidence antimicrobial resistance genes",
            'total_hits': len(df_filtered),
            'resistance_classes': resistance_classes,
            'top_hits': df_filtered.nlargest(10, 'bitscore')[['sseqid', 'pident', 'evalue', 'stitle']].to_dict('records')
        }
        
        # Save JSON report
        with open(output_dir / "antimicrobial_resistance_report.json", 'w') as f:
            import json
            json.dump(amr_report, f, indent=2)
        
        # Create visualization
        if resistance_classes:
            # Resistance class distribution
            class_counts = {k: len(v) for k, v in resistance_classes.items()}
            
            fig1 = px.bar(
                x=list(class_counts.keys()),
                y=list(class_counts.values()),
                title='Antimicrobial Resistance Gene Classes',
                labels={'x': 'Resistance Class', 'y': 'Number of Genes'},
                color=list(class_counts.values()),
                color_continuous_scale='Reds'
            )
            fig1.update_layout(template="plotly_white")
            fig1.write_html(output_dir / "amr_classes_distribution.html")
            
            # Identity distribution
            fig2 = px.histogram(
                df_filtered,
                x='pident',
                nbins=20,
                title='AMR Gene Identity Distribution',
                labels={'pident': 'Percentage Identity', 'count': 'Number of Genes'},
                color_discrete_sequence=['#FF6B6B']
            )
            fig2.update_layout(template="plotly_white")
            fig2.write_html(output_dir / "amr_identity_distribution.html")
        
        print(f"✓ Generated AMR report: {len(df_filtered)} resistance genes identified")
        
    except Exception as e:
        print(f"Error generating AMR report: {e}")


def generate_vf_report(vf_results, output_dir):
    """Generate virulence factor report from DIAMOND results"""
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    
    if not vf_results or not vf_results.exists():
        print("No virulence factor results to process")
        # Create empty report
        empty_report = {
            'summary': 'No virulence factors detected',
            'total_hits': 0,
            'virulence_categories': {},
            'top_hits': []
        }
        
        with open(output_dir / "virulence_factors_report.json", 'w') as f:
            import json
            json.dump(empty_report, f, indent=2)
        return
    
    try:
        # Read DIAMOND results
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle']
        df = pd.read_csv(vf_results, sep='\t', names=columns, header=None)
        
        if df.empty:
            print("No virulence factor hits found")
            return
        
        # Filter high-quality hits
        df_filtered = df[(df['pident'] >= 70) & (df['evalue'] <= 1e-10)]
        
        # Categorize virulence factors
        virulence_categories = {}
        vf_keywords = {
            'adhesion': ['adhesin', 'fimbri', 'pili', 'attach', 'binding'],
            'toxin': ['toxin', 'hemolysin', 'cytotoxin', 'enterotoxin'],
            'secretion_system': ['secretion', 'type III', 'type IV', 'type VI', 'T3SS', 'T4SS', 'T6SS'],
            'immune_evasion': ['capsule', 'LPS', 'immune', 'evasion', 'resistance'],
            'invasion': ['invasion', 'invasin', 'penetration'],
            'motility': ['flagell', 'motility', 'chemotaxis'],
            'regulation': ['regulator', 'sensor', 'response', 'quorum'],
            'stress_survival': ['stress', 'survival', 'persistence', 'dormancy']
        }
        
        for idx, row in df_filtered.iterrows():
            title_lower = row['stitle'].lower()
            categorized = False
            
            for category, keywords in vf_keywords.items():
                if any(keyword in title_lower for keyword in keywords):
                    if category not in virulence_categories:
                        virulence_categories[category] = []
                    virulence_categories[category].append({
                        'gene': row['sseqid'],
                        'identity': row['pident'],
                        'evalue': row['evalue'],
                        'description': row['stitle']
                    })
                    categorized = True
                    break
            
            if not categorized:
                if 'other' not in virulence_categories:
                    virulence_categories['other'] = []
                virulence_categories['other'].append({
                    'gene': row['sseqid'],
                    'identity': row['pident'],
                    'evalue': row['evalue'],
                    'description': row['stitle']
                })
        
        # Generate report
        vf_report = {
            'summary': f"Detected {len(df_filtered)} high-confidence virulence factors",
            'total_hits': len(df_filtered),
            'virulence_categories': virulence_categories,
            'top_hits': df_filtered.nlargest(10, 'bitscore')[['sseqid', 'pident', 'evalue', 'stitle']].to_dict('records')
        }
        
        # Save JSON report
        with open(output_dir / "virulence_factors_report.json", 'w') as f:
            import json
            json.dump(vf_report, f, indent=2)
        
        # Create visualization
        if virulence_categories:
            # Virulence category distribution
            category_counts = {k: len(v) for k, v in virulence_categories.items()}
            
            fig1 = px.pie(
                values=list(category_counts.values()),
                names=list(category_counts.keys()),
                title='Virulence Factor Categories',
                color_discrete_sequence=px.colors.qualitative.Set3
            )
            fig1.update_layout(template="plotly_white")
            fig1.write_html(output_dir / "virulence_categories_pie.html")
            
            # Identity vs E-value scatter plot
            fig2 = px.scatter(
                df_filtered,
                x='pident',
                y='evalue',
                size='bitscore',
                hover_data=['sseqid', 'stitle'],
                title='Virulence Factor Quality Metrics',
                labels={'pident': 'Percentage Identity', 'evalue': 'E-value'},
                log_y=True
            )
            fig2.update_layout(template="plotly_white")
            fig2.write_html(output_dir / "virulence_quality_scatter.html")
            
            # Top virulence factors table
            top_vf = df_filtered.nlargest(15, 'bitscore')
            fig3 = go.Figure(data=[go.Table(
                header=dict(values=['Gene ID', 'Identity %', 'E-value', 'Description'],
                           fill_color='lightblue',
                           align='left'),
                cells=dict(values=[
                    top_vf['sseqid'],
                    top_vf['pident'].round(1),
                    [f"{val:.2e}" for val in top_vf['evalue']],
                    [desc[:50] + '...' if len(desc) > 50 else desc for desc in top_vf['stitle']]
                ],
                fill_color='lavender',
                align='left'))
            ])
            fig3.update_layout(title="Top Virulence Factors")
            fig3.write_html(output_dir / "top_virulence_factors.html")
        
        print(f"✓ Generated virulence factor report: {len(df_filtered)} virulence genes identified")
        
    except Exception as e:
        print(f"Error generating virulence factor report: {e}")

def run_antimicrobial_resistance_scan(fasta_path, output_dir):
    """Scan for antimicrobial resistance genes using CARD database"""
    
    # Use CARD database files from your pathogen_db directory
    card_db = output_dir.parent / "databases" / "pathogen_db" / "protein_fasta_protein_homolog_model.fasta"
    
    if not card_db.exists():
        print("CARD database not found, skipping AMR scan")
        return None
    
    amr_out = output_dir / "amr_hits.txt"
    
    # Create temporary DIAMOND database if needed
    temp_dmnd = output_dir / "card_temp.dmnd"
    makedb_cmd = f"diamond makedb --in {card_db} --db {temp_dmnd.with_suffix('')}"
    
    try:
        subprocess.run(makedb_cmd, shell=True, check=True)
        
        # Run DIAMOND search
        cmd = f"diamond blastx -d {temp_dmnd.with_suffix('')} -q {fasta_path} -o {amr_out} " \
              "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle " \
              "--top 5 --evalue 1e-5 --threads 4"
        
        subprocess.run(cmd, shell=True, check=True)
        print("✓ AMR scan completed")
        
        # Clean up
        for ext in ['.dmnd']:
            temp_file = temp_dmnd.with_suffix(ext)
            if temp_file.exists():
                temp_file.unlink()
        
        return amr_out
        
    except subprocess.CalledProcessError as e:
        print(f"AMR scan failed: {e}")
        return None

def run_virulence_factor_scan(fasta_path, output_dir):
    """Scan for virulence factors using VFDB"""
    
    vfdb_file = output_dir.parent / "databases" / "pathogen_db" / "VFDB_setB_pro.fas"
    
    if not vfdb_file.exists():
        print("VFDB database not found, skipping virulence scan")
        return None
    
    vf_out = output_dir / "virulence_hits.txt"
    
    # Create temporary DIAMOND database
    temp_dmnd = output_dir / "vfdb_temp.dmnd"
    makedb_cmd = f"diamond makedb --in {vfdb_file} --db {temp_dmnd.with_suffix('')}"
    
    try:
        subprocess.run(makedb_cmd, shell=True, check=True)
        
        cmd = f"diamond blastx -d {temp_dmnd.with_suffix('')} -q {fasta_path} -o {vf_out} " \
              "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle " \
              "--top 5 --evalue 1e-5 --threads 4"
        
        subprocess.run(cmd, shell=True, check=True)
        print("✓ Virulence factor scan completed")
        
        # Clean up
        temp_dmnd.unlink()
        
        return vf_out
        
    except subprocess.CalledProcessError as e:
        print(f"Virulence scan failed: {e}")
        return None

def run_pathogen_scan(fasta_path, output_dir):
    """Screen for pathogens using Diamond BLAST against pathogen database"""
    
    # Define pathogen database path
    pathogen_db_path = Path("databases/pathogen_db/pathogen_db.dmnd")
    
    # Check if database exists
    if not pathogen_db_path.exists():
        print(f"Warning: Pathogen database not found at {pathogen_db_path}")
        print("Skipping pathogen screening...")
        # Create empty results file
        empty_results = output_dir / "pathogen_blast_results.txt"
        with open(empty_results, 'w') as f:
            f.write("# No pathogen database found\n")
        return empty_results
    
    blast_out = output_dir / "pathogen_blast_results.txt"
    
    print(f"Running pathogen screening against: {pathogen_db_path}")
    
    # Use simpler output format without taxonomy for now
    cmd = f"diamond blastx -d {pathogen_db_path.with_suffix('')} -q {fasta_path} -o {blast_out} " \
          "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle " \
          "--top 3 --evalue 1e-5 --threads 4"
    
    print(f"Running: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, timeout=1800)  # 30 min timeout
        print("✓ Pathogen screening completed")
        
        # Check if results file has content
        if blast_out.stat().st_size == 0:
            print("No pathogen hits found")
            with open(blast_out, 'w') as f:
                f.write("# No pathogen hits found\n")
        
        return blast_out
        
    except subprocess.TimeoutExpired:
        print("Warning: Pathogen screening timed out after 30 minutes")
        with open(blast_out, 'w') as f:
            f.write("# Pathogen screening timed out\n")
        return blast_out
        
    except subprocess.CalledProcessError as e:
        print(f"Warning: Pathogen screening failed: {e}")
        print("This might be due to database format issues")
        
        # Create empty results file so pipeline can continue
        with open(blast_out, 'w') as f:
            f.write("# Pathogen screening failed\n")
        return blast_out
    
    except Exception as e:
        print(f"Unexpected error in pathogen screening: {e}")
        with open(blast_out, 'w') as f:
            f.write("# Pathogen screening encountered an error\n")
        return blast_out