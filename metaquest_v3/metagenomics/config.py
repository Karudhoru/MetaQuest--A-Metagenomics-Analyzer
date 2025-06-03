from pathlib import Path
import plotly.express as px

# Configuration - Optimized for current database structure
DB_DIR = Path("./databases")

# KRAKEN2 TAXONOMIC DATABASE
KRAKEN_DB = DB_DIR  # Kraken files directly in databases/

# CAT DATABASE - NEW INTEGRATED PATHOGEN + TAXONOMY DATABASE
CAT_DB = DB_DIR / "cat_db" / "cat_database.dmnd"  # Primary pathogen database
CAT_TAXONOMY_MAP = DB_DIR / "cat_db" / "taxonomy.map"
CAT_FASTA_SOURCE = DB_DIR / "cat_db" / "pathogen_proteins.faa"
CAT_SEQID2TAXID = DB_DIR / "cat_db" / "seqid2taxid.map"
CAT_PATHOGEN_TAXIDS = DB_DIR / "cat_db" / "pathogen_taxids.txt"
CAT_PATHOGEN_SPECIES = DB_DIR / "cat_db" / "pathogen_species.txt"
CAT_AMR_CATEGORIES = DB_DIR / "cat_db" / "amr_categories.txt"

# LEGACY PATHOGEN DATABASES (for fallback)
PATHOGEN_DB = DB_DIR / "pathogen_db" / "pathogen_db_clean_tax.dmnd"
PATHOGEN_DB_WITH_TAX = DB_DIR / "pathogen_db" / "pathogen_db_with_tax.dmnd"

# SPECIALIZED DATABASES
SWISSPROT_DB = DB_DIR / "swissprot.dmnd"

# AMR AND VIRULENCE DATABASES
CARD_PROTEIN_DB = DB_DIR / "pathogen_db" / "protein_fasta_protein_homolog_model.fasta"
VFDB_DB = DB_DIR / "pathogen_db" / "VFDB_setB_pro.fas"
VFDB_PROCESSED = DB_DIR / "pathogen_db" / "vfdb_virulence.faa"

# TAXONOMY FILES (Multiple sources available)
TAXDUMP_DIR = DB_DIR / "taxdump"
TAXDUMP_CLEAN_DIR = DB_DIR / "taxdump_clean"  # Clean versions available
TAX_DIR = DB_DIR / "tax"  # Alternative tax directory

# Use the cleanest available taxonomy files
TAXDUMP_NODES = TAXDUMP_CLEAN_DIR / "nodes.dmp"
TAXDUMP_NAMES = TAXDUMP_CLEAN_DIR / "names.dmp"
PROT_ACCESSION2TAXID = TAXDUMP_DIR / "prot.accession2taxid.gz"

# PATHOGEN IDENTIFICATION - NOW PRE-COMPUTED
PATHOGEN_TAXIDS = CAT_PATHOGEN_TAXIDS  # Use CAT version
SEQID2TAXID_MAP = CAT_SEQID2TAXID     # Use CAT version

# VISUALIZATION COLORS
TAXONOMIC_COLORS = px.colors.qualitative.Set3
FUNCTIONAL_COLORS = px.colors.qualitative.Pastel
PATHOGEN_COLORS = {
    'Pathogenic': '#FF6B6B', 
    'Non-pathogenic': '#4ECDC4', 
    'Unknown': '#95A5A6',
    'AMR': '#E74C3C',
    'Virulence': '#F39C12'
}

    
BLAST_CONFIG = {
    'max_sequences_per_batch': 100,  # Limit to avoid API overload
    'rate_limit_delay': 0.4,         # 400ms between requests (NCBI recommendation)
    'max_retries': 3,                # Retry failed requests
    'timeout': 300,                  # 5 minute timeout per request
    'default_database': 'nt',        # Nucleotide database for FASTA
    'cache_expiry_days': 30          # Cache results for 30 days
}

# Cache Configuration
CACHE_CONFIG = {
    'enable_blast_cache': True,
    'cache_size_limit_mb': 500,      # Maximum cache size
    'cleanup_on_startup': True       # Clean old cache entries
}