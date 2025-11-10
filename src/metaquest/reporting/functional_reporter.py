"""
MetaQuest Functional Reporter
Enhanced functional annotation reporting with COG categories and transposase analysis
*** UPDATED TO USE GLOBAL OUTPUTFORMATTER ***
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter
import re
import json
from .base_reporter import BaseReporter
# (No need to import formatter, it's in BaseReporter)

class FunctionalReporter(BaseReporter):
    """
    Generate comprehensive functional annotation reports.
    Includes COG categories, feature composition, and annotation quality.
    """
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir)
        self.section_name = "Functional Report"
        
        # ... (COG category map) ...
        self.cog_categories = {
            'J': 'Translation, ribosomal structure',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'D': 'Cell cycle control, mitosis',
            'V': 'Defense mechanisms',
            'T': 'Signal transduction',
            'M': 'Cell wall/membrane biogenesis',
            'N': 'Cell motility',
            'U': 'Intracellular trafficking',
            'O': 'Post-translational modification, chaperones',
            'C': 'Energy production and conversion',
            'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis',
            'R': 'General function prediction only',
            'S': 'Function unknown',
            'X': 'Mobilome: transposases and prophages'
        }
    
    def generate_report(self,
                       sample_info_file: Path,
                       annotation_file: Path,
                       pathway_data_file: Optional[Path], # <-- Added this
                       functional_risk: Optional[Dict] = None) -> str:
        """
        Generate complete functional annotation report.
        """
        # --- Use formatter ---
        self.formatter.section_header("FUNCTIONAL ANNOTATION REPORT")

        stats = self._parse_sample_info(sample_info_file)
                
        # --- Use formatter ---
        self.formatter.debug(f"Func Report: Gene counts from sample.txt:")
        self.formatter.debug(f"   Total genes: {stats['gene']}")
        self.formatter.debug(f"   CDS: {stats['CDS']}")
        
        if stats['gene'] == 0:
            self.formatter.warning("Zero genes detected - functional report will be limited.")

        report_parts = []
        
        # 1. Genomic Features Overview
        report_parts.append(self._generate_features_overview(sample_info_file))

        if stats['gene'] == 0 or stats['CDS'] == 0:
            return '\n\n'.join(report_parts)
        
        # 2. Annotation Statistics
        report_parts.append(self._generate_annotation_stats(annotation_file, sample_info_file))
        
        # 3. COG Functional Categories
        report_parts.append(self._generate_cog_analysis(annotation_file))
        
        # 4. Transposase Analysis
        report_parts.append(self._generate_transposase_analysis(annotation_file))
        
        # 5. Top Annotated Functions
        report_parts.append(self._generate_top_functions(annotation_file))
        
        # 6. Annotation Quality Assessment
        report_parts.append(self._generate_quality_assessment(sample_info_file, annotation_file))
        
        # 7. Pathway Report
        report_parts.append(self._generate_pathway_report(pathway_data_file))
        
        return '\n\n'.join(report_parts)
    
    def _generate_features_overview(self, sample_file: Path) -> str:
        """Generate genomic features composition overview."""
        # Parse sample.txt
        stats = self._parse_sample_info(sample_file)
        
        total_genes = stats['gene']
        cds = stats['CDS']
        rrna = stats['rRNA']
        trna = stats['tRNA']
        
        lines = [
            "=" * 70,
            "🧬 GENOMIC FEATURES OVERVIEW",
            "=" * 70,
            "",
            "FEATURE COMPOSITION:",
            ""
        ]
        
        if total_genes == 0:
            lines.extend([
                "⚠️  No genes detected in assembly.",
                "",
                "POSSIBLE CAUSES:",
                "  • Contig filtering removed all sequences",
                "  • Minimum contig length too high",
                "  • Gene prediction failed",
                "  • Empty or corrupted assembly file",
                "  • sample.txt not updated properly",
                "",
                "TROUBLESHOOTING:",
                f"  1. Check Prokka directory: {sample_file.parent / 'prokka_annotation'}",
                "  2. Look for .tsv file with gene annotations",
                "  3. Verify assembly file has sequences",
                f"  4. Check if contigs were filtered (look for 'removed' in sample.txt)",
                "",
                "RECOMMENDATIONS:",
                "  1. Lower minimum contig length threshold (try --min-contig-length 200)",
                "  2. Disable contig filtering (--no-filter-contigs)",
                "  3. Check Prokka logs for errors",
                "  4. Verify input FASTA file integrity",
            ])
            return '\n'.join(lines) 
        
        # Calculate percentages (now safe because total_genes > 0)
        cds_pct = (cds / total_genes * 100)
        rrna_pct = (rrna / total_genes * 100)
        trna_pct = (trna / total_genes * 100)
        
        # Visual representation
        lines.extend([
            f"  CDS (Protein-coding):    {cds:>4} ({cds_pct:.1f}%)",
            f"  {'█' * int(cds_pct/2)}",
            f"",
            f"  rRNA (Ribosomal RNA):    {rrna:>4} ({rrna_pct:.1f}%)",
            f"  {'█' * int(rrna_pct/2)}",
            f"",
            f"  tRNA (Transfer RNA):     {trna:>4} ({trna_pct:.1f}%)",
            f"  {'█' * int(trna_pct/2)}",
            f"",
            f"  Total Genes:             {total_genes}",
        ])
        
        # Key metrics
        avg_gene_length = stats['bases'] // total_genes if total_genes > 0 else 0
        lines.extend([
            "",
            "KEY METRICS:",
            f"  Total Bases:       {stats['bases']:,} bp",
            f"  Number of Contigs: {stats['contigs']}",
            f"  Average Gene Length: {avg_gene_length:,} bp",
            f"  Coding Density:    {cds_pct:.1f}%",
        ])
        
        # Quality indicators
        lines.extend([
            "",
            "ASSEMBLY QUALITY INDICATORS:",
        ])
        
        # rRNA presence indicates completeness
        if rrna >= 3:
            lines.append("  ✅ rRNA genes present - suggests near-complete genome")
        else:
            lines.append("  ⚠️  Few rRNA genes - may be partial assembly")
        
        # tRNA presence
        if trna >= 15:
            lines.append("  ✅ Abundant tRNAs - good genome representation")
        elif trna > 0:
            lines.append("  ⚪ Some tRNAs present")
        else:
            lines.append("  ⚠️  No tRNAs - likely incomplete assembly")
        
        return '\n'.join(lines)

    
    def _generate_annotation_stats(self, annotation_file: Path, sample_file: Path) -> str:
        """Generate annotation coverage statistics from raw Diamond output."""
        lines = [
            "",
            "=" * 70,
            "📊 ANNOTATION STATISTICS",
            "=" * 70,
        ]
        
        # Load annotations (raw Diamond - no headers)
        if not annotation_file.exists() or annotation_file.stat().st_size == 0:
            return '\n'.join(lines + ["", "No functional annotations available"])
        
        try:
            df = pd.read_csv(annotation_file, sep='\t', header=None)
            
            # Assign standard Diamond column names
            diamond_cols = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                        'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                        'bitscore', 'description']
            df.columns = diamond_cols[:len(df.columns)]
        except Exception as e:
            return '\n'.join(lines + ["", f"Error loading annotations: {e}"])
        
        # Get sample stats
        stats = self._parse_sample_info(sample_file)
        total_cds = stats['CDS']
        
        # *** FIX: Check for zero CDS before division ***
        if total_cds == 0:
            lines.extend([
                "",
                "⚠️  No CDS (protein-coding sequences) found in assembly.",
                "   This may indicate:",
                "   • Assembly quality issues",
                "   • Filtering removed all contigs",
                "   • Gene prediction failed",
                "",
                "   Please check:",
                "   1. Prokka output logs",
                "   2. Contig filtering settings (minimum length)",
                "   3. Assembly statistics",
            ])
            return '\n'.join(lines)
        
        # Count annotated genes (unique query IDs)
        annotated_genes = df['query_id'].nunique()
        annotation_rate = (annotated_genes / total_cds * 100) if total_cds > 0 else 0
        
        # Count hypothetical proteins
        hypothetical = df[df['description'].str.contains('hypothetical', case=False, na=False)]
        hypothetical_count = hypothetical['query_id'].nunique()
        
        # *** FIX: Safe division with explicit check ***
        hypothetical_pct = (hypothetical_count / total_cds * 100) if total_cds > 0 else 0
        
        # SwissProt doesn't include COG IDs, so use functional categories
        functional_keywords = [
            'transposase', 'transport', 'permease', 'abc transporter',
            'metabolism', 'synthase', 'dehydrogenase', 'kinase',
            'replication', 'repair', 'transcription', 'polymerase',
            'ribosom', 'translation', 'trna', 'methyltransferase',
            'oxidoreductase', 'hydrolase', 'transferase', 'ligase',
            'secretion', 'membrane', 'binding'
        ]

        functionally_annotated = annotated_genes - hypothetical_count
        functional_pct = (functionally_annotated / total_cds * 100) if total_cds > 0 else 0

        lines.extend([
            "",
            f"Total CDS:              {total_cds}",
            f"Annotated CDS:          {annotated_genes} ({annotation_rate:.1f}%)",  # Uses full count
            f"Hypothetical proteins:  {hypothetical_count} ({hypothetical_pct:.1f}%)",
            f"Functionally annotated: {functionally_annotated} ({functional_pct:.1f}%)",  # Excludes hypothetical
            "",
            "ANNOTATION QUALITY:",
        ])

        # *** FIX: Update quality thresholds ***
        if annotation_rate >= 85:
            lines.append("  ✅ EXCELLENT coverage (≥85%)")
        elif annotation_rate >= 70:
            lines.append("  ✅ VERY GOOD coverage (70-85%)")
        elif annotation_rate >= 50:
            lines.append("  ✅ GOOD coverage (50-70%)")
        elif annotation_rate >= 30:
            lines.append("  ⚪ MODERATE coverage (30-50%)")
        else:
            lines.append("  ⚠️  LOW coverage (<30%)")
        
        # Average identity
        avg_identity = df['identity'].mean()
        lines.extend([
            "",
            f"Average sequence identity: {avg_identity:.1f}%",
        ])
        
        if avg_identity > 80:
            lines.append("  ✅ High similarity to known proteins")
        elif avg_identity > 50:
            lines.append("  ⚪ Moderate similarity")
        else:
            lines.append("  ⚠️  Low similarity - novel or divergent sequences")
        
        return '\n'.join(lines)
    
    def _generate_cog_analysis(self, annotation_file: Path) -> str:
        """Generate COG functional category analysis."""
        lines = [
            "",
            "=" * 70,
            "🔬 COG FUNCTIONAL CATEGORIES",
            "=" * 70,
        ]
        
        if not annotation_file.exists():
            return '\n'.join(lines + ["", "No annotations available"])
        
        df = pd.read_csv(annotation_file, sep='\t', header=None)
        df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                     'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                     'bitscore', 'description']
        
        # Extract COG IDs from descriptions
        cog_counts = Counter()
        category_counts = Counter()
        
        for desc in df['description']:
            # Look for COG pattern
            cog_matches = re.findall(r'COG\d{4}', str(desc))
            for cog in cog_matches:
                cog_counts[cog] += 1
            
            # Categorize by keywords (since we don't have explicit COG mapping)
            desc_lower = str(desc).lower()
            
            if any(word in desc_lower for word in ['transposase', 'transposon', 'insertion']):
                category_counts['X - Mobilome'] += 1
            elif any(word in desc_lower for word in ['ribosom', 'translation', 'trna']):
                category_counts['J - Translation'] += 1
            elif any(word in desc_lower for word in ['transport', 'permease', 'abc']):
                category_counts['P - Transport'] += 1
            elif any(word in desc_lower for word in ['metabolism', 'synthase', 'dehydrogenase']):
                category_counts['C/E/G - Metabolism'] += 1
            elif any(word in desc_lower for word in ['replication', 'dna', 'repair']):
                category_counts['L - Replication/Repair'] += 1
            elif any(word in desc_lower for word in ['transcription', 'rna polymerase']):
                category_counts['K - Transcription'] += 1
            elif 'hypothetical' in desc_lower:
                category_counts['S - Unknown function'] += 1
            else:
                category_counts['R - General function'] += 1
        
        # Display category breakdown
        total_categorized = sum(category_counts.values())
        
        lines.extend([
            "",
            "FUNCTIONAL CATEGORY DISTRIBUTION:",
            "",
        ])
        
        # Sort by count
        for category, count in sorted(category_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_categorized * 100) if total_categorized > 0 else 0
            bar = '█' * int(percentage / 2)
            lines.append(f"  {category:<30} {bar} {count:>3} ({percentage:>5.1f}%)")
        
        # Highlight key categories
        lines.extend([
            "",
            "KEY FUNCTIONAL HIGHLIGHTS:",
        ])
        
        if category_counts['X - Mobilome'] > 5:
            lines.append(
                f"  ⚠️  HIGH MOBILOME ACTIVITY: {category_counts['X - Mobilome']} mobile elements"
            )
            lines.append("      → Indicates active horizontal gene transfer potential")
        
        if category_counts['P - Transport'] > 10:
            lines.append(
                f"  ✅ ACTIVE TRANSPORT: {category_counts['P - Transport']} transport systems"
            )
            lines.append("      → Diverse nutrient acquisition capabilities")
        
        if category_counts['C/E/G - Metabolism'] > 15:
            lines.append(
                f"  ✅ METABOLIC DIVERSITY: {category_counts['C/E/G - Metabolism']} metabolic genes"
            )
        
        return '\n'.join(lines)
    
    def _generate_pathway_report(self, pathway_data_file: Optional[Path]) -> str:
        """
        Generates a human-readable metabolic pathway report from
        the pre-parsed EC count JSON file.
        """
        lines = [
            "",
            "=" * 70,
            "🔬 METABOLIC PATHWAY ANALYSIS (Lightweight)",
            "=" * 70,
        ]

        if not pathway_data_file or not pathway_data_file.exists():
            lines.append("\nNo EC number data found. Pathway analysis skipped.")
            return '\n'.join(lines)

        try:
            with open(pathway_data_file, 'r') as f:
                ec_counts = Counter(json.load(f))
        except Exception as e:
            lines.append(f"\nError loading pathway data: {e}")
            return '\n'.join(lines)

        if not ec_counts:
            lines.append("\nNo EC numbers found in annotations.")
            return '\n'.join(lines)

        lines.append(f"\nFound {len(ec_counts)} unique EC numbers with {sum(ec_counts.values())} total enzyme hits.")

        # This map stays in the reporter, as it's for presentation
        pathway_map = {
            'Glycolysis / Gluconeogenesis': ('1.1.1.1', '2.7.1.1', '5.3.1.9', '4.1.2.13', '1.2.1.12'),
            'TCA Cycle (Citrate Cycle)': ('4.1.3.7', '1.1.1.37', '1.2.4.1', '1.8.1.4', '2.3.1.12'),
            'Pentose Phosphate Pathway': ('1.1.1.49', '5.3.1.6', '3.1.1.31', '2.2.1.1', '5.1.3.1'),
            'Amino Acid Metabolism': ('1.4.1.', '2.6.1.', '4.1.1.', '6.1.1.'),
            'Fatty Acid Metabolism': ('1.1.1.35', '1.3.1.', '4.2.1.17', '1.1.1.211', '2.3.1.16'),
            'Nucleotide Metabolism': ('2.7.4.', '2.4.2.', '3.6.1.', '1.17.4.1'),
            'Antibiotic Biosynthesis': ('6.3.2.', '2.3.1.165', '1.1.1.238'),
        }

        pathway_counts = Counter()
        ec_to_pathway = {}

        for pathway, ec_prefixes in pathway_map.items():
            for ec, count in ec_counts.items():
                for prefix in ec_prefixes:
                    if ec.startswith(prefix):
                        pathway_counts[pathway] += count
                        ec_to_pathway[ec] = pathway
                        break

        lines.extend([
            "\n--- Pathway Summary ---",
            f"\n{'Metabolic Pathway':<40} | {'Gene Hits':>10}",
            f"{'-'*40:<40} | {'-'*10:>10}",
        ])

        for pathway, count in pathway_counts.most_common():
            lines.append(f"{pathway:<40} | {count:>10}")

        lines.extend([
            "\n\n--- Top 20 Most Abundant Enzymes ---",
            f"\n{'EC Number':<15} | {'Count':>7} | {'Pathway':<40}",
            f"{'-'*15:<15} | {'-'*7:>7} | {'-'*40:<40}",
        ])

        for ec, count in ec_counts.most_common(20):
            pathway = ec_to_pathway.get(ec, "Other")
            lines.append(f"{ec:<15} | {count:>7} | {pathway:<40}")

        return '\n'.join(lines)
    
    def _generate_transposase_analysis(self, annotation_file: Path) -> str:
        """Detailed analysis of transposases (critical for pathogenicity)."""
        lines = [
            "",
            "=" * 70,
            "🔄 MOBILE GENETIC ELEMENTS ANALYSIS",
            "=" * 70,
        ]
        
        if not annotation_file.exists():
            return '\n'.join(lines + ["", "No data available"])
        
        df = pd.read_csv(annotation_file, sep='\t', header=None)
        df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                     'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                     'bitscore', 'description']
        
        # Identify transposases
        transposase_keywords = ['transposase', 'transposon', 'IS element', 'IS66', 'IS30', 
                               'IS982', 'IS256', 'insertion sequence']
        
        transposases = df[df['description'].str.contains('|'.join(transposase_keywords), 
                                                        case=False, na=False)]
        
        if transposases.empty:
            lines.append("")
            lines.append("✅ No transposases detected")
            lines.append("   Low horizontal gene transfer potential")
            return '\n'.join(lines)
        
        transposon_count = len(transposases)
        unique_genes = transposases['query_id'].nunique()
        
        lines.extend([
            "",
            f"Total Transposase Hits: {transposon_count}",
            f"Unique Transposase Genes: {unique_genes}",
        ])
        
        # Classify IS families
        is_families = Counter()
        for desc in transposases['description']:
            if 'IS66' in desc:
                is_families['IS66'] += 1
            elif 'IS30' in desc:
                is_families['IS30'] += 1
            elif 'IS982' in desc:
                is_families['IS982'] += 1
            elif 'IS256' in desc:
                is_families['IS256'] += 1
            elif 'IS21' in desc:
                is_families['IS21'] += 1
            else:
                is_families['Other IS'] += 1
        
        lines.extend([
            "",
            "IS FAMILY DISTRIBUTION:",
        ])
        
        for family, count in sorted(is_families.items(), key=lambda x: x[1], reverse=True):
            lines.append(f"  ├─ {family:<15} {count:>3} genes")
        
        # Risk assessment
        lines.extend([
            "",
            "HORIZONTAL GENE TRANSFER RISK:",
        ])
        
        if unique_genes > 7:
            lines.extend([
                "  🔴 HIGH RISK",
                "     Multiple transposase families indicate active mobile",
                "     genetic elements. These can mobilize resistance and",
                "     virulence genes between bacteria.",
            ])
        elif unique_genes > 3:
            lines.extend([
                "  🟡 MODERATE RISK",
                "     Some mobile genetic elements present. Monitor for",
                "     acquisition of pathogenic traits.",
            ])
        else:
            lines.extend([
                "  🟢 LOW RISK",
                "     Few transposases detected.",
            ])
        
        return '\n'.join(lines)
    
    def _generate_top_functions(self, annotation_file: Path) -> str:
        """Display most abundant protein families."""
        lines = [
            "",
            "=" * 70,
            "🏆 TOP ANNOTATED FUNCTIONS",
            "=" * 70,
        ]
        
        if not annotation_file.exists():
            return '\n'.join(lines + ["", "No data available"])
        
        df = pd.read_csv(annotation_file, sep='\t', header=None)
        df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                     'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                     'bitscore', 'description']
        
        # Extract function names (simplified)
        def simplify_function(desc):
            # Remove organism names in brackets
            desc = re.sub(r'\[.*?\]', '', desc)
            # Remove database IDs
            desc = re.sub(r'[A-Z]{2,}_\d+', '', desc)
            # Take first part before detailed info
            return desc.split('[')[0].strip()[:50]
        
        df['function'] = df['description'].apply(simplify_function)
        
        # Count occurrences
        function_counts = df['function'].value_counts().head(15)
        
        lines.extend([
            "",
            f"{'Function':<52} {'Count':>8}",
            "-" * 65,
        ])
        
        for func, count in function_counts.items():
            lines.append(f"{func[:50]:<52} {count:>8}")
        
        return '\n'.join(lines)
    
    def _generate_quality_assessment(self, sample_file: Path, annotation_file: Path) -> str:
        """Overall annotation quality assessment."""
        lines = [
            "",
            "=" * 70,
            "✅ ANNOTATION QUALITY ASSESSMENT",
            "=" * 70,
        ]
        
        stats = self._parse_sample_info(sample_file)
        
        if not annotation_file.exists():
            lines.append("")
            lines.append("⚠️  No functional annotations available")
            return '\n'.join(lines)
        
        df = pd.read_csv(annotation_file, sep='\t', header=None)
        df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                     'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                     'bitscore', 'description']
        
        annotated = df['query_id'].nunique()
        total_cds = stats['CDS']
        annotation_rate = (annotated / total_cds * 100) if total_cds > 0 else 0
        
        # Calculate quality score (0-100)
        quality_score = 0
        
        # Annotation coverage (40 points)
        if annotation_rate > 70:
            quality_score += 40
        elif annotation_rate > 50:
            quality_score += 30
        elif annotation_rate > 30:
            quality_score += 20
        else:
            quality_score += 10
        
        # Average identity (30 points)
        avg_identity = df['identity'].mean()
        if avg_identity > 80:
            quality_score += 30
        elif avg_identity > 60:
            quality_score += 20
        elif avg_identity > 40:
            quality_score += 10
        
        # Functional diversity (30 points)
        unique_functions = df['description'].nunique()
        if unique_functions > 30:
            quality_score += 30
        elif unique_functions > 20:
            quality_score += 20
        elif unique_functions > 10:
            quality_score += 10
        
        lines.extend([
            "",
            f"OVERALL QUALITY SCORE: {quality_score}/100",
            "",
            "BREAKDOWN:",
            f"  Annotation Coverage:   {annotation_rate:.1f}% {'✅' if annotation_rate > 50 else '⚠️'}",
            f"  Average Identity:      {avg_identity:.1f}% {'✅' if avg_identity > 60 else '⚠️'}",
            f"  Functional Diversity:  {unique_functions} unique functions {'✅' if unique_functions > 20 else '⚠️'}",
            "",
            "INTERPRETATION:",
        ])
        
        if quality_score >= 80:
            lines.append("  🌟 EXCELLENT - High-quality comprehensive annotations")
        elif quality_score >= 60:
            lines.append("  ✅ GOOD - Reliable functional annotation")
        elif quality_score >= 40:
            lines.append("  ⚪ MODERATE - Acceptable but could be improved")
        else:
            lines.append("  ⚠️  LOW - Consider additional annotation methods")
        
        return '\n'.join(lines)
    
    def _parse_sample_info(self, sample_file: Path) -> Dict:
        """Parse sample.txt file for assembly statistics."""
        stats = {
            'organism': '', 'contigs': 0, 'bases': 0, 'CDS': 0,
            'gene': 0, 'rRNA': 0, 'tRNA': 0, 'removed_contigs': 0
        }
        
        if not sample_file.exists():
            # --- Use formatter ---
            self.formatter.warning(f"Sample info file not found: {sample_file}")
            return stats
        
        with open(sample_file) as f:
            for line in f:
                line = line.strip()
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # FIX: Use case-insensitive matching
                    key_lower = key.lower()
                    
                    if key_lower == 'organism':
                        stats['organism'] = value
                    elif key_lower == 'cds':
                        stats['CDS'] = int(value)
                    elif key_lower == 'gene':
                        stats['gene'] = int(value)
                    elif key_lower == 'rrna':
                        stats['rRNA'] = int(value)
                    elif key_lower == 'trna':
                        stats['tRNA'] = int(value)
                    elif key_lower == 'contigs':
                        stats['contigs'] = int(value)
                    elif key_lower == 'bases':
                        stats['bases'] = int(value)
                    elif 'removed' in key_lower or 'filtered' in key_lower:
                        stats['removed_contigs'] = int(value)
        
        # FIX: If gene count is 0 but we have a Prokka directory, try reading from there
        if stats['gene'] == 0 and stats['CDS'] == 0:
            # --- Use formatter ---
            self.formatter.warning("No genes in sample.txt - attempting fallback")
            prokka_stats = self._load_prokka_stats_fallback(sample_file.parent)
            if prokka_stats:
                stats.update(prokka_stats)
        
        return stats
    
    def _load_prokka_stats_fallback(self, output_dir: Path) -> Dict:
        """
        Fallback: Read gene counts directly from Prokka TSV file.
        """
        stats = {}
        prokka_dir = output_dir / "prokka_annotation"
        
        if not prokka_dir.exists():
            # --- Use formatter ---
            self.formatter.warning(f"Prokka directory not found: {prokka_dir}")
            return stats
        
        tsv_files = list(prokka_dir.glob("*.tsv"))
        if not tsv_files:
            # --- Use formatter ---
            self.formatter.warning(f"No TSV files found in {prokka_dir}")
            return stats
        
        tsv_file = tsv_files[0]
        # --- Use formatter ---
        self.formatter.info(f"Reading Prokka annotations from: {tsv_file.name}")
        
        try:
            # Read Prokka TSV (tab-separated, with header)
            df = pd.read_csv(tsv_file, sep='\t')
            
            # Count feature types
            if 'ftype' in df.columns:
                feature_counts = df['ftype'].value_counts()
                
                stats['CDS'] = int(feature_counts.get('CDS', 0))
                stats['rRNA'] = int(feature_counts.get('rRNA', 0))
                stats['tRNA'] = int(feature_counts.get('tRNA', 0))
                
                # Count unique locus tags for total genes
                if 'locus_tag' in df.columns:
                    stats['gene'] = df['locus_tag'].nunique()
                else:
                    stats['gene'] = len(df[df['ftype'] == 'gene'])
                
                self.formatter.success(f"Loaded from Prokka TSV: {stats['gene']} genes, {stats['CDS']} CDS")
            else:
                self.formatter.error(f"Failed to read Prokka TSV: {e}")
            
            # *** ADD THIS: Get base count from assembly file ***
            # Look for the FASTA file used for annotation
            fasta_files = list(prokka_dir.glob("*.fna"))  # Prokka's nucleotide output
            if not fasta_files:
                # Try filtered assembly in parent directory
                fasta_files = list(output_dir.glob("*_filtered.fasta"))
            
            if fasta_files:
                from Bio import SeqIO
                total_bases = sum(len(record.seq) for record in SeqIO.parse(fasta_files[0], 'fasta'))
                stats['bases'] = total_bases
                self.formatter.info(f"✓ Calculated bases from assembly: {total_bases:,} bp")
            
            # Get contig count from filter log
            filter_log = output_dir / "contig_filtering.log"
            if filter_log.exists():
                with open(filter_log, 'r') as f:
                    for line in f:
                        if 'Filtered contigs:' in line:
                            stats['contigs'] = int(line.split(':')[1].strip())
                        
        except Exception as e:
            self.formatter.error(f"Failed to read Prokka TSV: {e}")
        
        return stats
    
    def generate_visualization_data(self,
                                   sample_file: Path,
                                   annotation_file: Path) -> Dict:
        """
        Generate data structure for visualization module.
        
        Returns:
            Dict with data ready for functional_visualizer.py
        """
        stats = self._parse_sample_info(sample_file)
        
        viz_data = {
            'features': {
                'CDS': stats['CDS'],
                'rRNA': stats['rRNA'],
                'tRNA': stats['tRNA'],
                'total_genes': stats['gene']
            },
            'assembly': {
                'total_bases': stats['bases'],
                'contigs': stats['contigs'],
                'avg_gene_length': stats['bases'] // stats['gene'] if stats['gene'] > 0 else 0
            }
        }
        
        if annotation_file.exists():
            df = pd.read_csv(annotation_file, sep='\t', header=None)
            df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                         'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                         'bitscore', 'description']
            
            # COG categories
            category_counts = self._extract_cog_categories(df)
            viz_data['cog_categories'] = category_counts
            
            # Transposases
            transposases = df[df['description'].str.contains('transposase', case=False, na=False)]
            viz_data['transposases'] = {
                'count': len(transposases),
                'unique_genes': transposases['query_id'].nunique()
            }
            
            # Annotation quality
            viz_data['annotation_quality'] = {
                'annotated_count': df['query_id'].nunique(),
                'annotation_rate': df['query_id'].nunique() / stats['CDS'] if stats['CDS'] > 0 else 0,
                'avg_identity': float(df['identity'].mean())
            }
        
        return viz_data
    
    def _extract_cog_categories(self, df: pd.DataFrame) -> Dict:
        """Extract COG category counts from annotations."""
        category_counts = Counter()
        
        for desc in df['description']:
            desc_lower = str(desc).lower()
            
            if any(word in desc_lower for word in ['transposase', 'transposon']):
                category_counts['Mobilome'] += 1
            elif any(word in desc_lower for word in ['ribosom', 'translation']):
                category_counts['Translation'] += 1
            elif any(word in desc_lower for word in ['transport', 'permease', 'abc']):
                category_counts['Transport'] += 1
            elif any(word in desc_lower for word in ['metabolism', 'synthase']):
                category_counts['Metabolism'] += 1
            elif any(word in desc_lower for word in ['replication', 'dna']):
                category_counts['Replication'] += 1
            elif any(word in desc_lower for word in ['transcription', 'rna polymerase']):
                category_counts['Transcription'] += 1
            elif 'hypothetical' in desc_lower:
                category_counts['Unknown'] += 1
            else:
                category_counts['Other'] += 1
        
        return dict(category_counts)
    
    def export_annotations_table(self, annotation_file: Path, output_file: Path):
        """Export detailed annotations table to CSV."""
        if not annotation_file.exists():
            self.formatter.warning(f"Annotation file not found: {annotation_file}")
            return
        
        df = pd.read_csv(annotation_file, sep='\t', header=None)
        df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                     'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                     'bitscore', 'description']
        
        # Add simplified function column
        df['function'] = df['description'].apply(lambda x: re.sub(r'\[.*?\]', '', str(x)).strip()[:100])
        
        # Add category column
        df['category'] = df['description'].apply(self._categorize_function)
        
        # Export
        df.to_csv(output_file, index=False)
        self.formatter.success(f"Annotations table exported to {output_file.name}")
    
    def _categorize_function(self, description: str) -> str:
        """Categorize function into broad category."""
        desc_lower = str(description).lower()
        
        if any(word in desc_lower for word in ['transposase', 'transposon']):
            return 'Mobilome'
        elif any(word in desc_lower for word in ['ribosom', 'translation']):
            return 'Translation'
        elif any(word in desc_lower for word in ['transport', 'permease', 'abc']):
            return 'Transport'
        elif any(word in desc_lower for word in ['metabolism', 'synthase', 'dehydrogenase']):
            return 'Metabolism'
        elif any(word in desc_lower for word in ['replication', 'dna', 'repair']):
            return 'Replication/Repair'
        elif any(word in desc_lower for word in ['transcription', 'rna polymerase']):
            return 'Transcription'
        elif 'hypothetical' in desc_lower:
            return 'Unknown Function'
        else:
            return 'Other'