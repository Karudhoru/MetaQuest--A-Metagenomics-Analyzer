# metagenomics/file_validator.py

import os
from pathlib import Path
from Bio import SeqIO
import gzip
import hashlib
import numpy as np
from collections import Counter
from datetime import datetime

class FileValidator:
    """Comprehensive file validation and statistics for MetaQuest"""
    
    def __init__(self):
        self.min_sequences = 100
        self.quality_threshold = 20
        self.overrep_threshold = 10.0  # 10% threshold for stopping pipeline
        self.adapter_threshold = 5.0   # 5% adapter content threshold
        self.q30_threshold = 80.0      # 80% Q30 bases minimum
        self.sample_size = 10000       # Sample first N sequences for contamination check
        self.common_adapters = [
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "CTGTCTCTTATACACATCT",
        ]
        
    def validate_and_analyze(self, input_file, file_type):
        """Main validation and analysis function with improved status reporting."""
        print(f"\n{'─'*70}")
        print(f"📋 FILE VALIDATION REPORT")
        print(f"{'─'*70}")
        print(f"File:     {Path(input_file).name}")
        print(f"Path:     {input_file}")
        print(f"Type:     {file_type.upper()}")
        print(f"Time:     {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'─'*70}\n")
        
        if not self._check_file_exists(input_file): 
            return False, None
        if not self._validate_format(input_file, file_type): 
            return False, None
            
        print("⏳ Analyzing file contents...\n")
        stats = self._generate_statistics(input_file, file_type)
        
        self._display_statistics(stats, file_type)
        
        has_errors, error_messages = self._get_validation_status(stats, file_type)
        
        print(f"\n{'─'*70}")
        if has_errors:
            print("❌ VALIDATION STATUS: FAILED")
            print(f"{'─'*70}")
            print("\n🛑 CRITICAL ISSUES DETECTED:")
            for msg in error_messages:
                print(f"   • {msg}")
            print("\n💡 RECOMMENDATIONS:")
            self._print_recommendations(stats, file_type, error_messages)
            print(f"\n{'─'*70}\n")
            return False, stats
        else:
            print("✅ VALIDATION STATUS: PASSED")
            print(f"{'─'*70}")
            print("\n✓ All quality criteria met - file is ready for analysis")
            print(f"{'─'*70}\n")
            return True, stats

    def _generate_statistics(self, input_file, file_type):
        """Generate comprehensive file statistics"""
        stats = {
            'file_path': input_file,
            'file_size_mb': os.path.getsize(input_file) / (1024 * 1024),
            'is_compressed': input_file.endswith('.gz'),
            'md5_checksum': self._calculate_md5(input_file)[:16]
        }
        
        if file_type == 'fastq':
            stats.update(self._analyze_fastq(input_file))
        else:
            stats.update(self._analyze_fasta(input_file))
            
        return stats

    def _analyze_fastq(self, fastq_file):
        """Optimized FASTQ analysis with single-pass algorithm."""
        stats = {
            'total_sequences': 0, 'total_bases': 0, 'min_length': float('inf'),
            'max_length': 0, 'mean_length': 0, 'median_length': 0,
            'n50_length': 0, 'gc_content': 0, 'mean_quality': 0,
            'min_quality': float('inf'), 'max_quality': 0, 'q20_bases': 0,
            'q30_bases': 0, 'quality_encoding': None,
            'adapter_content_percent': 0, 'overrepresented_sequences': [],
            'top_5_sequences': []
        }
        
        handle = gzip.open(fastq_file, 'rt') if fastq_file.endswith('.gz') else open(fastq_file, 'r')
        
        lengths = []
        gc_count = 0
        total_quality_sum = 0
        sequence_counts = Counter()
        adapter_hits = 0
        sampled_sequences = 0
        
        try:
            print("   [Stage 1/2] Reading sequences...", end='', flush=True)
            for i, record in enumerate(SeqIO.parse(handle, 'fastq')):
                stats['total_sequences'] += 1
                seq_str = str(record.seq)
                seq_len = len(seq_str)
                
                lengths.append(seq_len)
                stats['min_length'] = min(stats['min_length'], seq_len)
                stats['max_length'] = max(stats['max_length'], seq_len)
                
                stats['total_bases'] += seq_len
                gc_count += seq_str.count('G') + seq_str.count('C')
                
                qual_scores = record.letter_annotations.get('phred_quality', [])
                if qual_scores:
                    total_quality_sum += sum(qual_scores)
                    stats['min_quality'] = min(stats['min_quality'], *qual_scores)
                    stats['max_quality'] = max(stats['max_quality'], *qual_scores)
                    stats['q20_bases'] += sum(1 for q in qual_scores if q >= 20)
                    stats['q30_bases'] += sum(1 for q in qual_scores if q >= 30)

                # Sample first N sequences for contamination check
                if sampled_sequences < self.sample_size:
                    sequence_counts[seq_str] += 1
                    for adapter in self.common_adapters:
                        if adapter in seq_str:
                            adapter_hits += 1
                            break
                    sampled_sequences += 1
            
            print(" Done!")
            print("   [Stage 2/2] Computing statistics...", end='', flush=True)
        finally:
            handle.close()

        # Calculate derived statistics
        if lengths:
            stats['mean_length'] = np.mean(lengths)
            stats['median_length'] = np.median(lengths)
            stats['n50_length'] = self._calculate_n50(lengths)
        
        if stats['total_bases'] > 0:
            stats['gc_content'] = (gc_count / stats['total_bases']) * 100
            stats['mean_quality'] = total_quality_sum / stats['total_bases']
            
            if stats['min_quality'] < 59 and stats['max_quality'] <= 74:
                stats['quality_encoding'] = 'Phred+33 (Sanger/Illumina 1.8+)'
            else:
                stats['quality_encoding'] = 'Phred+64 (Illumina 1.3-1.5)'

        if sampled_sequences > 0:
            stats['adapter_content_percent'] = (adapter_hits / sampled_sequences) * 100
            for seq, count in sequence_counts.most_common(5):
                percentage = (count / sampled_sequences) * 100
                stats['top_5_sequences'].append((seq, count, percentage))
                if percentage >= self.overrep_threshold:
                    stats['overrepresented_sequences'].append((seq, count, percentage))

        print(" Done!\n")
        return stats
    
    def _analyze_fasta(self, fasta_file):
        """Optimized FASTA analysis."""
        stats = {
            'total_sequences': 0, 'total_bases': 0, 'min_length': float('inf'),
            'max_length': 0, 'mean_length': 0, 'median_length': 0,
            'gc_content': 0, 'ambiguous_bases': 0, 'duplicate_ids': 0,
            'sequence_type': 'unknown', 'avg_gc_per_sequence': 0,
            'gc_variance': 0, 'length_variance': 0, 'n_count': 0, 'gap_count': 0,
        }
        
        handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file, 'r')
        
        lengths = []
        gc_percentages = []
        gc_count = 0
        seen_ids = set()
        protein_chars = set('DEFHIKLMNPQRSVWY')
        nucleotide_only = set('ATCGUN')
        
        try:
            print("   [Stage 1/2] Reading sequences...", end='', flush=True)
            for i, record in enumerate(SeqIO.parse(handle, 'fasta')):
                seq_str = str(record.seq).upper()
                seq_len = len(seq_str)
                
                if record.id in seen_ids:
                    stats['duplicate_ids'] += 1
                seen_ids.add(record.id)
                
                stats['total_sequences'] += 1
                stats['total_bases'] += seq_len
                lengths.append(seq_len)
                
                stats['min_length'] = min(stats['min_length'], seq_len)
                stats['max_length'] = max(stats['max_length'], seq_len)
                
                seq_gc = seq_str.count('G') + seq_str.count('C')
                gc_count += seq_gc
                if seq_len > 0:
                    gc_percentages.append((seq_gc / seq_len) * 100)
                
                stats['ambiguous_bases'] += sum(1 for base in seq_str if base not in 'ATCG')
                stats['n_count'] += seq_str.count('N')
                stats['gap_count'] += seq_str.count('-')
                
                if i < 10 and stats['sequence_type'] == 'unknown':
                    unique_chars = set(seq_str.replace('-', '').replace('*', ''))
                    if unique_chars.intersection(protein_chars):
                        stats['sequence_type'] = 'protein'
                    elif unique_chars.issubset(nucleotide_only):
                        stats['sequence_type'] = 'nucleotide'
            
            print(" Done!")
            print("   [Stage 2/2] Computing statistics...", end='', flush=True)
        finally:
            handle.close()
        
        if lengths:
            stats['mean_length'] = np.mean(lengths)
            stats['median_length'] = np.median(lengths)
            stats['length_variance'] = np.var(lengths)
            if stats['mean_length'] > 0:
                stats['length_cv'] = (np.std(lengths) / stats['mean_length']) * 100
        
        if stats['total_bases'] > 0 and stats['sequence_type'] != 'protein':
            stats['gc_content'] = (gc_count / stats['total_bases']) * 100
            
        if gc_percentages and stats['sequence_type'] != 'protein':
            stats['avg_gc_per_sequence'] = np.mean(gc_percentages)
            stats['gc_variance'] = np.var(gc_percentages)
            stats['gc_std'] = np.std(gc_percentages)
        
        if stats['sequence_type'] == 'protein':
            stats['sequence_category'] = 'Protein sequences'
        elif stats['total_sequences'] == 1:
            stats['sequence_category'] = 'Single sequence (genome/chromosome)'
        elif stats['mean_length'] > 10000:
            stats['sequence_category'] = 'Assembled contigs/scaffolds'
        elif stats['mean_length'] > 1000:
            stats['sequence_category'] = 'Gene sequences'
        else:
            stats['sequence_category'] = 'Short sequences'
        
        print(" Done!\n")
        return stats
    
    def _display_statistics(self, stats, file_type):
        """Display statistics in a professional scientific format."""
        
        # File Information
        print("📊 SUMMARY STATISTICS")
        print("═" * 70)
        print(f"\n{'File Properties':<30} {'Value':<40}")
        print("─" * 70)
        print(f"{'  Size':<30} {stats['file_size_mb']:.2f} MB")
        print(f"{'  Compression':<30} {'gzip' if stats['is_compressed'] else 'None'}")
        print(f"{'  MD5 Checksum':<30} {stats['md5_checksum']}...")
        
        if file_type == 'fastq':
            self._display_fastq_stats(stats)
        else:
            self._display_fasta_stats(stats)
    
    def _display_fastq_stats(self, stats):
        """Display FASTQ statistics in scientific format."""
        
        print(f"\n{'Sequence Metrics':<30} {'Value':<40}")
        print("─" * 70)
        print(f"{'  Total Sequences':<30} {stats['total_sequences']:,}")
        print(f"{'  Total Bases':<30} {stats['total_bases']:,} bp")
        print(f"{'  Length (min/mean/max)':<30} {stats['min_length']}/{stats['mean_length']:.0f}/{stats['max_length']} bp")
        print(f"{'  Median Length':<30} {stats['median_length']:.0f} bp")
        print(f"{'  N50':<30} {stats['n50_length']:.0f} bp")
        print(f"{'  GC Content':<30} {stats['gc_content']:.1f}%")
        
        print(f"\n{'Quality Metrics':<30} {'Value':<40}")
        print("─" * 70)
        print(f"{'  Encoding':<30} {stats['quality_encoding'] or 'Unknown'}")
        print(f"{'  Mean Quality Score':<30} {stats['mean_quality']:.1f}")
        print(f"{'  Quality Range':<30} Q{stats['min_quality']:.0f} - Q{stats['max_quality']:.0f}")
        
        q20_pct = (stats['q20_bases']/stats['total_bases']*100)
        q30_pct = (stats['q30_bases']/stats['total_bases']*100)
        print(f"{'  Bases ≥ Q20':<30} {q20_pct:.1f}%")
        print(f"{'  Bases ≥ Q30':<30} {q30_pct:.1f}%")
        
        print(f"\n{'Contamination Analysis':<30} {'Value':<40}")
        print("─" * 70)
        print(f"{'  Adapter Content':<30} {stats['adapter_content_percent']:.2f}% (sampled)")
        
        if stats['overrepresented_sequences']:
            print(f"{'  Overrepresented Seqs':<30} {len(stats['overrepresented_sequences'])} detected (≥{self.overrep_threshold}%)")
            print(f"\n  {'Top Overrepresented Sequences':<55} {'Frequency':<15}")
            print(f"  {'-'*55:<55} {'-'*15:<15}")
            for seq, count, pct in stats['overrepresented_sequences'][:5]:
                seq_display = seq[:52] + '...' if len(seq) > 52 else seq
                print(f"  {seq_display:<55} {pct:>6.2f}%")
        else:
            print(f"{'  Overrepresented Seqs':<30} None detected")
    
    def _display_fasta_stats(self, stats):
        """Display FASTA statistics in scientific format."""
        
        print(f"\n{'Sequence Metrics':<30} {'Value':<40}")
        print("─" * 70)
        print(f"{'  Sequence Category':<30} {stats['sequence_category']}")
        print(f"{'  Total Sequences':<30} {stats['total_sequences']:,}")
        print(f"{'  Total Bases':<30} {stats['total_bases']:,}")
        print(f"{'  Length (min/mean/max)':<30} {stats['min_length']}/{stats['mean_length']:.0f}/{stats['max_length']}")
        print(f"{'  Median Length':<30} {stats['median_length']:.0f}")
        
        if stats['sequence_type'] != 'protein':
            print(f"\n{'Composition Analysis':<30} {'Value':<40}")
            print("─" * 70)
            print(f"{'  GC Content':<30} {stats['gc_content']:.1f}%")
            if stats['total_sequences'] > 1:
                print(f"{'  GC per Sequence (μ±σ)':<30} {stats['avg_gc_per_sequence']:.1f}% ± {stats.get('gc_std', 0):.1f}%")
            if stats['n_count'] > 0:
                n_pct = (stats['n_count'] / stats['total_bases']) * 100
                print(f"{'  Ambiguous Bases (N)':<30} {stats['n_count']:,} ({n_pct:.2f}%)")
            if stats['gap_count'] > 0:
                gap_pct = (stats['gap_count'] / stats['total_bases']) * 100
                print(f"{'  Gap Characters':<30} {stats['gap_count']:,} ({gap_pct:.2f}%)")
        
        if stats['total_sequences'] > 1:
            print(f"\n{'Sequence Diversity':<30} {'Value':<40}")
            print("─" * 70)
            print(f"{'  Length CV':<30} {stats.get('length_cv', 0):.1f}%")
            if stats['duplicate_ids'] > 0:
                print(f"{'  Duplicate IDs':<30} {stats['duplicate_ids']} detected")
    
    def _get_validation_status(self, stats, file_type):
        """Check for critical errors that should stop the pipeline."""
        errors = []

        if file_type == 'fastq':
            # Critical thresholds
            if stats['total_sequences'] < self.min_sequences:
                errors.append(f"Insufficient sequences: {stats['total_sequences']} < {self.min_sequences} required")
            
            if stats['mean_quality'] < self.quality_threshold:
                errors.append(f"Low quality score: Q{stats['mean_quality']:.1f} < Q{self.quality_threshold} required")
            
            q30_pct = (stats.get('q30_bases', 0) / stats.get('total_bases', 1) * 100)
            if q30_pct < self.q30_threshold:
                errors.append(f"Low Q30 content: {q30_pct:.1f}% < {self.q30_threshold}% required")
            
            if stats.get('adapter_content_percent', 0) > self.adapter_threshold:
                errors.append(f"High adapter contamination: {stats['adapter_content_percent']:.2f}% > {self.adapter_threshold}% threshold")
            
            if stats.get('overrepresented_sequences'):
                max_overrep = max(pct for _, _, pct in stats['overrepresented_sequences'])
                errors.append(f"Overrepresented sequences detected: {len(stats['overrepresented_sequences'])} seq (max {max_overrep:.1f}%)")
        
        else:  # fasta
            if stats['total_sequences'] == 0:
                errors.append("No valid sequences found in file")
            
            if stats['duplicate_ids'] > 0:
                errors.append(f"Duplicate sequence IDs: {stats['duplicate_ids']} duplicates found")

        return len(errors) > 0, errors
    
    def _print_recommendations(self, stats, file_type, error_messages):
        """Print actionable recommendations based on errors."""
        
        if file_type == 'fastq':
            print("\n   Quality preprocessing is required before analysis.")
            print("   Recommended tools:")
            print()
            print("   1. fastp (recommended)")
            print("      fastp -i input.fq -o output.fq --qualified_quality_phred 20 \\")
            print("           --detect_adapter_for_pe --cut_right")
            print()
            print("   2. Trimmomatic")
            print("      trimmomatic SE input.fq output.fq ILLUMINACLIP:adapters.fa:2:30:10 \\")
            print("           LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50")
            print()
            print("   3. cutadapt")
            print("      cutadapt -a ADAPTER -q 20 -m 50 -o output.fq input.fq")
            print()
            print("   After preprocessing, rerun validation to verify improvements.")
        else:
            print("\n   Please address the following issues:")
            if stats.get('duplicate_ids', 0) > 0:
                print("   • Remove or rename duplicate sequence IDs")
            if stats.get('total_sequences', 0) == 0:
                print("   • Verify file format and content")
    
    def _calculate_n50(self, lengths):
        """Calculate N50 metric."""
        if not lengths:
            return 0
        sorted_lengths = sorted(lengths, reverse=True)
        cumsum = np.cumsum(sorted_lengths)
        total = cumsum[-1]
        for i, cs in enumerate(cumsum):
            if cs >= total / 2:
                return sorted_lengths[i]
        return 0
    
    def _calculate_md5(self, filepath):
        """Calculate MD5 checksum efficiently."""
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def _check_file_exists(self, filepath):
        """Check if file exists and is readable."""
        if not Path(filepath).exists():
            print(f"\n❌ ERROR: File not found")
            print(f"   Path: {filepath}")
            return False
        if not os.access(filepath, os.R_OK):
            print(f"\n❌ ERROR: File not readable")
            print(f"   Path: {filepath}")
            return False
        return True
    
    def _validate_format(self, filepath, expected_type):
        """Validate file format quickly."""
        try:
            handle = gzip.open(filepath, 'rt') if filepath.endswith('.gz') else open(filepath, 'r')
            
            if expected_type == 'fastq':
                next(SeqIO.parse(handle, 'fastq'))
            else:
                next(SeqIO.parse(handle, 'fasta'))
            
            handle.close()
            return True
            
        except StopIteration:
            print(f"\n❌ ERROR: File is empty or contains no valid sequences")
            return False
        except Exception as e:
            print(f"\n❌ ERROR: Invalid {expected_type.upper()} format")
            print(f"   Details: {str(e)}")
            return False