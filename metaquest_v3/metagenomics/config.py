from pathlib import Path
import plotly.express as px

# Configuration - Updated paths based on your directory structure
DB_DIR = Path("./databases")  # Relative to script location
KRAKEN_DB = DB_DIR  # The kraken files are directly in databases/
SWISSPROT_DB = DB_DIR/"swissprot.dmnd"
# Primary pathogen database (currently lacking taxonomy - needs rebuilding)
PATHOGEN_DB = DB_DIR / "pathogen_db" / "pathogen_db.dmnd"  # Current database without taxonomy

# NEW: Pathogen database with taxonomy (to be created)
PATHOGEN_DB_WITH_TAX = DB_DIR / "pathogen_db" / "pathogen_db_with_tax.dmnd"  # Fixed version

# Source FASTA file for rebuilding the database
PATHOGEN_FASTA_SOURCE = DB_DIR / "pathogen_db" / "protein_fasta_protein_homolog_model.fasta"

# Alternative pathogen databases
PATHOGEN_MARKERS_DB = DB_DIR / "pathogen_db" / "pathogen_markers.dmnd"  # Alternative pathogen markers
PATHOGEN_FASTA_DB = DB_DIR / "pathogen_db" / "pathogen_markers.faa"  # FASTA version

# TAXONOMY CONFIGURATION - REQUIRED FOR PATHOGEN SCANNING
TAXDUMP_DIR = DB_DIR / "taxdump"
# These files are needed for Diamond taxonomy integration:
TAXDUMP_NODES = TAXDUMP_DIR / "nodes.dmp"
TAXDUMP_NAMES = TAXDUMP_DIR / "names.dmp" 
PROT_ACCESSION2TAXID = TAXDUMP_DIR / "prot.accession2taxid.gz"  # From NCBI

# PATHOGEN IDENTIFICATION FILES
PATHOGEN_TAXIDS = DB_DIR / "pathogenic_taxids.txt"  # You'll need to create this
SEQID2TAXID_MAP = DB_DIR / "seqid2taxid.map"  # For taxonomy mapping

# AMR AND VIRULENCE DATABASES (Referenced in the document)
CARD_DB = DB_DIR / "pathogen_db" / "protein_fasta_protein_homolog_model.fasta"  # CARD AMR database
VFDB_DB = DB_DIR / "pathogen_db" / "VFDB_setB_pro.fas"  # Virulence factor database

# VISUALIZATION COLORS
TAXONOMIC_COLORS = px.colors.qualitative.Set3
FUNCTIONAL_COLORS = px.colors.qualitative.Pastel
PATHOGEN_COLORS = {'Pathogenic': '#FF6B6B', 'Non-pathogenic': '#4ECDC4', 'Unknown': '#95A5A6'}