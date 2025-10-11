#!/bin/bash
# MetaQuest Advanced Database Setup Script
#
# This script downloads and formats all necessary databases for the MetaQuest tool.
# It includes checks for existing files and supports selective downloads.

set -e # Exit immediately if a command exits with a non-zero status.

# --- Configuration ---
DB_DIR="$(dirname "$0")/../databases"

# --- Helper Functions ---
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "❌ ERROR: Command not found: '$1'. Please install it first."
        echo "   (Hint: 'conda install -c bioconda $1')"
        exit 1
    fi
}

# --- Main Logic ---
main() {
    echo "--- MetaQuest Database Setup ---"
    
    # 1. Check for required command-line tools
    echo "1. Checking for required tools..."
    check_command wget
    check_command tar
    check_command diamond
    echo "✅ All tools found."

    # 2. Parse command-line arguments
    if [ "$#" -eq 0 ]; then
        echo "No specific database requested. Use --all to get everything, or see --help."
        exit 0
    fi
    
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --kraken) setup_kraken=true; shift ;;
            --pathogen) setup_pathogen=true; shift ;;
            --swissprot) setup_swissprot=true; shift ;;
            --all) setup_kraken=true; setup_pathogen=true; setup_swissprot=true; shift ;;
            --help) show_help; exit 0 ;;
            *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
        esac
    done

    # 3. Create base directory
    mkdir -p "$DB_DIR"
    echo -e "\nDatabase directory: $DB_DIR"

    # 4. Run setup functions based on flags
    if [ "$setup_kraken" = true ]; then
        setup_kraken_db
    fi
    if [ "$setup_pathogen" = true ]; then
        setup_pathogen_dbs
    fi
    if [ "$setup_swissprot" = true ]; then
        setup_swissprot_db
    fi

    echo -e "\n🎉 Database setup complete!"
}

show_help() {
    echo "Usage: ./setup_databases.sh [options]"
    echo
    echo "Options:"
    echo "  --kraken       Download and set up the Kraken2/Bracken database."
    echo "  --pathogen     Download and set up the CARD and VFDB pathogen databases."
    echo "  --swissprot    Download and set up the SwissProt+COG combined database for functional annotation."
    echo "  --all          Download and set up all available databases."
    echo "  --help         Show this help message."
}

setup_kraken_db() {
    echo -e "\n--- Setting up Kraken2/Bracken Database ---"
    if [ -f "$DB_DIR/hash.k2d" ]; then
        echo "✅ Kraken2 database already exists. Skipping."
    else
        echo "Downloading MiniKraken database (~8GB)..."
        wget -q --show-progress -nc -P "$DB_DIR" https://genome-idx.s3.amazonaws.com/kraken/minikraken_8GB_202003.tgz
        echo "Decompressing..."
        tar -xzf "$DB_DIR/minikraken_8GB_202003.tgz" -C "$DB_DIR" --strip-components=1
        rm "$DB_DIR/minikraken_8GB_202003.tgz"
        echo "✅ Kraken2 setup complete."
    fi
}

setup_pathogen_dbs() {
    echo -e "\n--- Setting up Pathogen Databases (CARD & VFDB) ---"
    PATHOGEN_DIR="$DB_DIR/pathogen_db"
    mkdir -p "$PATHOGEN_DIR"

    # CARD
    if [ ! -f "$PATHOGEN_DIR/protein_fasta_protein_homolog_model.fasta" ]; then
        echo "Downloading CARD database..."
        wget -q --show-progress -nc -P "$PATHOGEN_DIR" https://card.mcmaster.ca/latest/data
        tar -xjf "$PATHOGEN_DIR/data" -C "$PATHOGEN_DIR" protein_fasta_protein_homolog_model.fasta
        rm "$PATHOGEN_DIR/data"
    else
        echo "✅ CARD database file already exists. Skipping download."
    fi

    # VFDB
    if [ ! -f "$PATHOGEN_DIR/VFDB_setB_pro.fas" ]; then
        echo "Downloading VFDB..."
        wget -q --show-progress -nc -P "$PATHOGEN_DIR" http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas
    else
        echo "✅ VFDB file already exists. Skipping download."
    fi
    echo "✅ Pathogen databases downloaded."
}

setup_swissprot_db() {
    echo -e "\n--- Setting up SwissProt + COG Combined Database ---"
    
    # Final combined database name matching config.py
    COMBINED_DB="$DB_DIR/SwissProt_COG_db.dmnd"
    
    if [ -f "$COMBINED_DB" ]; then
        echo "✅ SwissProt+COG combined DIAMOND database already exists. Skipping."
    else
        # Download SwissProt
        if [ ! -f "$DB_DIR/uniprot_sprot.fasta" ]; then
            echo "Downloading SwissProt FASTA file..."
            wget -q --show-progress -nc -P "$DB_DIR" https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz
            echo "Decompressing SwissProt..."
            gunzip "$DB_DIR/uniprot_sprot.fasta.gz"
        else
            echo "✅ SwissProt FASTA already exists."
        fi
        
        # Download COG
        if [ ! -f "$DB_DIR/cog.fasta" ]; then
            echo "Downloading COG database..."
            wget -q --show-progress -nc -P "$DB_DIR" https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
            echo "Decompressing COG..."
            gunzip "$DB_DIR/cog-20.fa.gz"
            mv "$DB_DIR/cog-20.fa" "$DB_DIR/cog.fasta"
        else
            echo "✅ COG FASTA already exists."
        fi
        
        # Combine both databases
        echo "Combining SwissProt and COG databases..."
        cat "$DB_DIR/uniprot_sprot.fasta" "$DB_DIR/cog.fasta" > "$DB_DIR/combined_swissprot_cog.fasta"
        
        # Format combined database with DIAMOND
        echo "Formatting combined database with DIAMOND (this may take several minutes)..."
        diamond makedb --in "$DB_DIR/combined_swissprot_cog.fasta" -d "$DB_DIR/SwissProt_COG_db"
        
        # Cleanup intermediate files
        echo "Cleaning up intermediate files..."
        rm "$DB_DIR/uniprot_sprot.fasta" "$DB_DIR/cog.fasta" "$DB_DIR/combined_swissprot_cog.fasta"
        
        echo "✅ SwissProt+COG combined database setup complete."
    fi
}

# Run the main function with all command line arguments
main "$@"