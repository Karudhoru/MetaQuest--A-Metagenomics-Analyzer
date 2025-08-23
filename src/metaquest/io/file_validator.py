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
        self.min_length = 50
        self.quality_threshold = 20
        self.overrep_threshold = 0.1  # <-- NEW: Default threshold
        self.common_adapters = [
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "CTGTCTCTTATACACATCT",
        ]
        
    def validate_and_analyze(self, input_file, file_type):
        """Main validation and analysis function"""
        print(f"\n{'='*60}")
        print(f"🔍 METAQUEST FILE VALIDATION")
        print(f"{'='*60}")
        print(f"File: {input_file}")
        print(f"Type: {file_type.upper()}")
        print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*60}\n")
        
        # Check file exists and is readable
        if not self._check_file_exists(input_file):
            return False, None
            
        # Check file format
        if not self._validate_format(input_file, file_type):
            return False, None
            
        # Generate comprehensive statistics
        print("⏳ Analyzing file contents...")
        stats = self._generate_statistics(input_file, file_type)
        
        # Display statistics
        self._display_statistics(stats, file_type)
        
        # --- MODIFIED STATUS LOGIC ---
        has_errors, has_warnings = self._get_validation_status(stats, file_type)
        
        print(f"\n{'='*60}")
        if has_errors:
            print("❌ FILE VALIDATION: FAILED")
            print("   Please address the critical issues noted in the 'VALIDATION CRITERIA' section.")
        elif has_warnings:
            print("⚠️ FILE VALIDATION: PASSED WITH WARNINGS")
            print("   The file is technically valid, but please review the recommendations before proceeding with analysis.")
        else: # No errors, no warnings
            print("✅ FILE VALIDATION: PASSED")
        print(f"{'='*60}\n")
        
        return not has_errors, stats
    
    def _generate_statistics(self, input_file, file_type):
        """Generate comprehensive file statistics"""
        stats = {
            'file_path': input_file,
            'file_size_mb': os.path.getsize(input_file) / (1024 * 1024),
            'file_type': file_type,
            'is_compressed': input_file.endswith('.gz'),
            'md5_checksum': self._calculate_md5(input_file)[:16]  # First 16 chars
        }
        
        if file_type == 'fastq':
            stats.update(self._analyze_fastq(input_file))
        else:  # fasta
            stats.update(self._analyze_fasta(input_file))
            
        return stats
    
    def _get_validation_status(self, stats, file_type):
        """Checks for both critical errors and non-critical warnings."""
        has_errors = False
        has_warnings = False

        # --- 1. Check for critical errors that should stop the analysis ---
        if file_type == 'fastq':
            if stats['total_sequences'] < self.min_sequences:
                has_errors = True
            if stats['mean_quality'] < self.quality_threshold:
                has_errors = True
        
        # --- 2. Check for non-critical warnings ---
        if stats.get('overrepresented_sequences'):
            has_warnings = True
        if stats.get('adapter_content_percent', 0) > 5.0: # Flag a warning for >5% adapter content
            has_warnings = True

        return has_errors, has_warnings
    
    def _analyze_fastq(self, fastq_file):
        """Detailed FASTQ analysis with a correct and efficient single-pass algorithm."""
        stats = {
            'total_sequences': 0, 'total_bases': 0, 'min_length': float('inf'),
            'max_length': 0, 'mean_length': 0, 'median_length': 0,
            'n50_length': 0, 'gc_content': 0, 'mean_quality': 0,
            'min_quality': float('inf'), 'max_quality': 0, 'q20_bases': 0,
            'q30_bases': 0, 'ambiguous_bases': 0,
            'adapter_content_percent': 0, 'overrepresented_sequences': [],
            'quality_encoding': None,
            'top_5_sequences': [],
        }
        
        handle = gzip.open(fastq_file, 'rt') if fastq_file.endswith('.gz') else open(fastq_file, 'r')
        
        lengths, qualities = [], []
        gc_count = 0
        sequence_counts = Counter()
        adapter_hits = 0
        
        try:
            print("  Processing: ", end='', flush=True)
            # --- NEW SIMPLIFIED SINGLE-PASS LOOP ---
            # This loop iterates through the entire file and performs all checks.
            for i, record in enumerate(SeqIO.parse(handle, 'fastq')):
                if i > 0 and i % 20000 == 0:
                    print(".", end='', flush=True)
                
                stats['total_sequences'] += 1
                seq_str = str(record.seq)
                seq_len = len(seq_str)
                
                # Full-file stats
                lengths.append(seq_len)
                stats['min_length'] = min(stats['min_length'], seq_len)
                stats['max_length'] = max(stats['max_length'], seq_len)
                
                qual_scores = record.letter_annotations.get('phred_quality', [])
                if qual_scores:
                    qualities.extend(qual_scores)
                    stats['q20_bases'] += sum(1 for q in qual_scores if q >= 20)
                    stats['q30_bases'] += sum(1 for q in qual_scores if q >= 30)

                stats['total_bases'] += seq_len
                gc_count += seq_str.count('G') + seq_str.count('C')
                
                # Duplication and Adapter stats
                sequence_counts[seq_str] += 1
                for adapter in self.common_adapters:
                    if adapter in seq_str:
                        adapter_hits += 1
                        break
            
            print(" Done!")
        finally:
            handle.close()

        # --- Calculate all derived statistics AFTER the loop ---
        if lengths:
            stats['mean_length'] = np.mean(lengths)
            stats['median_length'] = np.median(lengths)
            stats['n50_length'] = self._calculate_n50(lengths)
        if stats['total_bases'] > 0:
            stats['gc_content'] = (gc_count / stats['total_bases']) * 100
        
        if qualities:
            stats['mean_quality'] = np.mean(qualities)
            stats['min_quality'] = min(qualities)
            stats['max_quality'] = max(qualities)
            if min(qualities) < 59 and max(qualities) <= 74:
                stats['quality_encoding'] = 'Sanger/Illumina 1.8+ (Phred+33)'
            else:
                stats['quality_encoding'] = 'Illumina 1.3-1.5 (Phred+64)'

        # Calculate duplication and adapter stats based on the TOTAL sequences
        if stats['total_sequences'] > 0:
            stats['adapter_content_percent'] = (adapter_hits / stats['total_sequences']) * 100
            
            for seq, count in sequence_counts.most_common(15):
                percentage = (count / stats['total_sequences']) * 100
                if len(stats['top_5_sequences']) < 5:
                    stats['top_5_sequences'].append((seq, count, percentage))
                
                if percentage >= self.overrep_threshold:
                    stats['overrepresented_sequences'].append((seq, count, percentage))

        return stats
    
    def _analyze_fasta(self, fasta_file):
        """Detailed FASTA analysis - focused on sequence characteristics"""
        stats = {
            'total_sequences': 0,
            'total_bases': 0,
            'min_length': float('inf'),
            'max_length': 0,
            'mean_length': 0,
            'median_length': 0,
            'gc_content': 0,
            'ambiguous_bases': 0,
            'duplicate_ids': 0,
            'sequence_type': 'unknown',  # protein/nucleotide
            'avg_gc_per_sequence': 0,
            'gc_variance': 0,
            'length_variance': 0,
            'n_count': 0,  # Count of N bases
            'gap_count': 0,  # Count of gaps (-)
        }
        
        handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file, 'r')
        
        lengths = []
        gc_percentages = []
        gc_count = 0
        seen_ids = set()
        protein_chars = set('DEFHIKLMNPQRSVWY')
        nucleotide_only = set('ATCGUN')
        
        try:
            # Process sequences with progress indicator
            print("  Processing: ", end='', flush=True)
            for i, record in enumerate(SeqIO.parse(handle, 'fasta')):
                if i % 1000 == 0:
                    print(".", end='', flush=True)
                
                seq_str = str(record.seq).upper()
                seq_len = len(seq_str)
                
                # Check for duplicate IDs
                if record.id in seen_ids:
                    stats['duplicate_ids'] += 1
                seen_ids.add(record.id)
                
                # Update basic stats
                stats['total_sequences'] += 1
                stats['total_bases'] += seq_len
                lengths.append(seq_len)
                
                # Length stats
                stats['min_length'] = min(stats['min_length'], seq_len)
                stats['max_length'] = max(stats['max_length'], seq_len)
                
                # GC content (for nucleotides)
                seq_gc = seq_str.count('G') + seq_str.count('C')
                gc_count += seq_gc
                if seq_len > 0:
                    gc_percentages.append((seq_gc / seq_len) * 100)
                
                # Check sequence content
                stats['ambiguous_bases'] += sum(1 for base in seq_str if base not in 'ATCG')
                stats['n_count'] += seq_str.count('N')
                stats['gap_count'] += seq_str.count('-')
                
                # Detect sequence type (first 10 sequences)
                if i < 10 and stats['sequence_type'] == 'unknown':
                    unique_chars = set(seq_str.replace('-', '').replace('*', ''))
                    if unique_chars.intersection(protein_chars):
                        stats['sequence_type'] = 'protein'
                    elif unique_chars.issubset(nucleotide_only):
                        stats['sequence_type'] = 'nucleotide'
            
            print(" Done!")
        
        finally:
            handle.close()
        
        # Calculate derived statistics
        if lengths:
            stats['mean_length'] = np.mean(lengths)
            stats['median_length'] = np.median(lengths)
            stats['length_variance'] = np.var(lengths)
            
            # Coefficient of variation for length
            if stats['mean_length'] > 0:
                stats['length_cv'] = (np.std(lengths) / stats['mean_length']) * 100
        
        if stats['total_bases'] > 0 and stats['sequence_type'] != 'protein':
            stats['gc_content'] = (gc_count / stats['total_bases']) * 100
            
        if gc_percentages and stats['sequence_type'] != 'protein':
            stats['avg_gc_per_sequence'] = np.mean(gc_percentages)
            stats['gc_variance'] = np.var(gc_percentages)
            stats['gc_std'] = np.std(gc_percentages)
        
        # Determine likely sequence category
        if stats['sequence_type'] == 'protein':
            stats['sequence_category'] = 'Protein sequences'
        elif stats['total_sequences'] == 1:
            stats['sequence_category'] = 'Single sequence (likely genome/chromosome)'
        elif stats['mean_length'] > 10000:
            stats['sequence_category'] = 'Assembled contigs/scaffolds'
        elif stats['mean_length'] > 1000:
            stats['sequence_category'] = 'Gene sequences'
        else:
            stats['sequence_category'] = 'Short sequences'
        
        return stats
    
    def _display_statistics(self, stats, file_type):
            """Display statistics in a beautiful shell format"""
            # File Information
            print("\n📁 FILE INFORMATION")
            print("─" * 40)
            print(f"  Size:        {stats['file_size_mb']:.2f} MB")
            print(f"  MD5:         {stats['md5_checksum']}...")
            print(f"  Compressed:  {'Yes (gzip)' if stats['is_compressed'] else 'No'}")
            
            if file_type == 'fastq':
                print("\n🔬 CONTAMINATION & DUPLICATION")
                print("─" * 40)
                print(f"  Adapter Content:    {stats['adapter_content_percent']:.2f}% of reads sampled")
                
                # Display the flagged sequences in a clean table
                if stats['overrepresented_sequences']:
                    print(f"\n  Flagged as Overrepresented (>{self.overrep_threshold:.2f}%) - Top 10:")
                    print(f"  {'Sequence':<40} {'Percentage'}")
                    print(f"  {'-'*40:<40} {'-'*10}")
                    # Loop through the top 10 flagged sequences
                    for seq, count, percentage in stats['overrepresented_sequences'][:10]:
                        print(f"  {seq[:37]+'...':<40} {percentage:.4f}%")
                else:
                    print("\n  Overrepresented Seq:  None found above threshold.")
                # --- END OF UPDATED SECTION ---

                # RECOMMENDATIONS
                if stats['adapter_content_percent'] > 5 or stats['overrepresented_sequences']:
                    print("\n  RECOMMENDATIONS:")
                    if stats['adapter_content_percent'] > 5:
                        print("  - ⚠️ High adapter content detected. Pre-processing with a tool like 'fastp' or 'cutadapt' is strongly recommended.")
                    if stats['overrepresented_sequences']:
                        print("  - ℹ️ Overrepresented sequences found. This may indicate PCR duplication. Consider investigating or using a deduplication tool.")
                
                self._display_fastq_stats(stats)
            else:
                self._display_fasta_stats(stats)
    
    def _display_fastq_stats(self, stats):
        """Display FASTQ-specific statistics."""
        # Sequence Statistics
        print("\n🧬 SEQUENCE STATISTICS")
        print("─" * 40)
        print(f"  Total Sequences:  {stats['total_sequences']:,}")
        print(f"  Total Bases:      {stats['total_bases']:,}")
        print(f"  Length Range:     {stats['min_length']:,} - {stats['max_length']:,} bp")
        print(f"  Mean Length:      {stats['mean_length']:.1f} bp")
        print(f"  Median Length:    {stats['median_length']:.1f} bp")
        print(f"  N50 Length:       {stats['n50_length']:.1f} bp")
        print(f"  GC Content:       {stats['gc_content']:.1f}%")
        
        # Quality Statistics
        print("\n📊 QUALITY STATISTICS")
        print("─" * 40)
        print(f"  Encoding:         {stats['quality_encoding'] or 'Unknown'}")
        print(f"  Mean Quality:     {stats['mean_quality']:.1f}")
        print(f"  Quality Range:    {stats['min_quality']} - {stats['max_quality']}")
        print(f"  Q20 Bases:        {(stats['q20_bases']/stats['total_bases']*100):.1f}%")
        print(f"  Q30 Bases:        {(stats['q30_bases']/stats['total_bases']*100):.1f}%")
        
        # Quality Assessment
        print("\n🔍 QUALITY ASSESSMENT")
        print("─" * 40)

        warnings = []

        if stats.get('min_quality') == stats.get('max_quality'):
            print(f"  Note: All bases have uniform quality ({stats.get('min_quality')})")
            
        if stats['mean_quality'] < 20:
            warnings.append(("Low mean quality", f"{stats['mean_quality']:.1f} < 20", "⚠️"))
        
        if stats.get('ambiguous_bases', 0) > stats.get('total_bases', 1) * 0.01:
            pct = (stats['ambiguous_bases']/stats['total_bases']*100)
            warnings.append(("High ambiguous bases", f"{pct:.2f}%", "⚠️"))
        
        # --- THIS IS THE FIX ---
        # Check against the new percentage metric instead of the old boolean key
        if stats.get('adapter_content_percent', 0) > 1.0:
            warnings.append((f"High Adapter Content ({stats['adapter_content_percent']:.2f}%)", "Detected", "⚠️"))
        
        if (stats.get('q30_bases', 0) / stats.get('total_bases', 1) * 100) < 80:
            pct = (stats['q30_bases']/stats['total_bases']*100)
            warnings.append(("Low Q30 percentage", f"{pct:.1f}% < 80%", "⚠️"))
        
        if warnings:
            for issue, value, icon in warnings:
                print(f"  {icon} {issue}: {value}")
        else:
            print("  ✅ All quality metrics passed")
    
    def _display_fasta_stats(self, stats):
        """Display FASTA-specific statistics - more relevant metrics"""
        # Sequence Statistics
        print("\n🧬 SEQUENCE STATISTICS")
        print("─" * 40)
        print(f"  Type:             {stats['sequence_category']}")
        print(f"  Total Sequences:  {stats['total_sequences']:,}")
        if stats['sequence_type'] == 'protein':
            print(f"  Total Residues:   {stats['total_bases']:,}")
        else:
            print(f"  Total Bases:      {stats['total_bases']:,}")
        print(f"  Length Range:     {stats['min_length']:,} - {stats['max_length']:,}")
        print(f"  Mean Length:      {stats['mean_length']:.1f}")
        print(f"  Median Length:    {stats['median_length']:.1f}")
        
        # Only show GC content for nucleotide sequences
        if stats['sequence_type'] != 'protein':
            print("\n🧪 COMPOSITION ANALYSIS")
            print("─" * 40)
            print(f"  GC Content:       {stats['gc_content']:.1f}%")
            if stats['total_sequences'] > 1:
                print(f"  GC per Sequence:  {stats['avg_gc_per_sequence']:.1f}% ± {stats.get('gc_std', 0):.1f}%")
            if stats['n_count'] > 0:
                n_pct = (stats['n_count'] / stats['total_bases']) * 100
                print(f"  N Bases:          {stats['n_count']:,} ({n_pct:.2f}%)")
            if stats['gap_count'] > 0:
                gap_pct = (stats['gap_count'] / stats['total_bases']) * 100
                print(f"  Gaps:             {stats['gap_count']:,} ({gap_pct:.2f}%)")
        
        # Sequence Diversity (for multi-sequence files)
        if stats['total_sequences'] > 1:
            print("\n📊 SEQUENCE DIVERSITY")
            print("─" * 40)
            print(f"  Length CV:        {stats.get('length_cv', 0):.1f}%")
            if stats['sequence_type'] != 'protein' and stats.get('gc_variance', 0) > 0:
                print(f"  GC Variation:     {stats.get('gc_std', 0):.1f}%")
        
        # Quality Assessment
        print("\n🔍 QUALITY ASSESSMENT")
        print("─" * 40)
        
        warnings = []
        if stats['duplicate_ids'] > 0:
            warnings.append(("Duplicate IDs found", f"{stats['duplicate_ids']} duplicates", "⚠️"))
        
        if stats['sequence_type'] != 'protein' and stats['ambiguous_bases'] > stats['total_bases'] * 0.01:
            pct = (stats['ambiguous_bases']/stats['total_bases']*100)
            warnings.append(("High ambiguous bases", f"{pct:.2f}%", "⚠️"))
        
        if stats['total_sequences'] == 1 and stats['mean_length'] < 1000:
            warnings.append(("Very short sequence", f"{stats['mean_length']:.0f} bp", "⚠️"))
        
        if stats['total_sequences'] > 1000 and stats.get('length_cv', 0) > 200:
            warnings.append(("Highly variable lengths", f"CV: {stats.get('length_cv', 0):.1f}%", "⚠️"))
        
        if warnings:
            for issue, value, icon in warnings:
                print(f"  {icon} {issue}: {value}")
        else:
            print("  ✅ All quality checks passed")
    
    def _calculate_n50(self, lengths):
        """Calculate N50 metric"""
        if not lengths:
            return 0
        sorted_lengths = sorted(lengths, reverse=True)
        cumsum = np.cumsum(sorted_lengths)
        total = cumsum[-1]
        for i, cs in enumerate(cumsum):
            if cs >= total / 2:
                return sorted_lengths[i]
        return 0
    
    def _calculate_n90(self, lengths):
        """Calculate N90 metric"""
        if not lengths:
            return 0
        sorted_lengths = sorted(lengths, reverse=True)
        cumsum = np.cumsum(sorted_lengths)
        total = cumsum[-1]
        for i, cs in enumerate(cumsum):
            if cs >= total * 0.9:
                return sorted_lengths[i]
        return 0
    
    def _calculate_md5(self, filepath):
        """Calculate MD5 checksum of file"""
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            # Read in chunks to handle large files
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def _check_file_exists(self, filepath):
        """Check if file exists and is readable"""
        if not Path(filepath).exists():
            print(f"\n❌ ERROR: File not found: {filepath}")
            return False
        if not os.access(filepath, os.R_OK):
            print(f"\n❌ ERROR: File not readable: {filepath}")
            return False
        return True
    
    def _validate_format(self, filepath, expected_type):
        """Validate file format"""
        try:
            handle = gzip.open(filepath, 'rt') if filepath.endswith('.gz') else open(filepath, 'r')
            
            # Try to parse first record
            if expected_type == 'fastq':
                next(SeqIO.parse(handle, 'fastq'))
            else:
                next(SeqIO.parse(handle, 'fasta'))
            
            handle.close()
            return True
            
        except StopIteration:
            print(f"\n❌ ERROR: File is empty or has no valid sequences")
            return False
        except Exception as e:
            print(f"\n❌ ERROR: Invalid {expected_type.upper()} format")
            print(f"   Details: {str(e)}")
            return False
    
    def _validate_quality(self, stats, file_type):
        """Validate file quality based on statistics"""
        print("\n⚙️  VALIDATION CRITERIA")
        print("─" * 40)
        
        is_valid = True
        
        if file_type == 'fastq':
            # Check sequence count
            print(f"  Minimum sequences:  {self.min_sequences:,}", end='')
            if stats['total_sequences'] < self.min_sequences:
                print(f" ❌ (Found: {stats['total_sequences']:,})")
                is_valid = False
            else:
                print(f" ✅ (Found: {stats['total_sequences']:,})")
            
            # Check quality
            print(f"  Minimum quality:    {self.quality_threshold}", end='')
            if stats['mean_quality'] < self.quality_threshold:
                print(f" ❌ (Found: {stats['mean_quality']:.1f})")
                is_valid = False
            else:
                print(f" ✅ (Found: {stats['mean_quality']:.1f})")
                
        else:  # fasta
            # Check sequence count
            print(f"  Minimum sequences:  1", end='')
            if stats['total_sequences'] == 0:
                print(f" ❌ (Found: 0)")
                is_valid = False
            else:
                print(f" ✅ (Found: {stats['total_sequences']:,})")
            
            # Check for duplicates
            print(f"  Unique IDs:         Required", end='')
            if stats['duplicate_ids'] > 0:
                print(f" ⚠️  (Duplicates will be renamed)")
            else:
                print(f" ✅ (All unique)")
        
        return is_valid