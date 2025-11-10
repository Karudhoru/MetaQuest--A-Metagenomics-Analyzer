#!/bin/bash
# MetaQuest Database Setup Script (Updated)
#
# This script downloads MiniKraken, builds a custom pathogen DB,
# and builds a combined SwissProt+COG DB.

set -e # Exit immediately if a command exits with a non-zero status.

# --- Configuration ---
DB_DIR="$(dirname "$0")/../databases"

# --- ⬇️ USER: CONFIGURE YOUR CUSTOM PATHOGEN DB HERE ---
#
# Name of your custom FASTA file (must be inside $DB_DIR)
CUSTOM_PATHOGEN_FASTA_NAME="metaquest_pathogen_markers.faa"
#
# Desired name for the output DIAMOND database (will create "pathogen_markers.dmnd")
CUSTOM_PATHOGEN_DB_NAME="metaquest_pathogen_markers"
#
# --- End of Configuration ---


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
    check_command gunzip
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
        setup_custom_pathogen_db
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
    echo "  --kraken       Download and set up the Kraken2 (MiniKraken) database."
    echo "  --pathogen     Build a DIAMOND database from your custom pathogen FASTA."
    echo "                 (Requires $CUSTOM_PATHOGEN_FASTA_NAME to be in $DB_DIR)"
    echo "  --swissprot    Download and set up the SwissProt+COG combined database."
    echo "  --all          Run all setup steps (Kraken, Pathogen, SwissProt)."
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

setup_custom_pathogen_db() {
    echo -e "\n--- Setting up Custom Pathogen Database ---"
    
    FASTA_IN="$DB_DIR/$CUSTOM_PATHOGEN_FASTA_NAME"
    DB_OUT_BASE="$DB_DIR/$CUSTOM_PATHOGEN_DB_NAME" # diamond adds .dmnd
    DB_OUT_FINAL="$DB_OUT_BASE.dmnd"

    # Check if final DB already exists
    if [ -f "$DB_OUT_FINAL" ]; then
        echo "✅ Custom pathogen DIAMOND database already exists: $DB_OUT_FINAL. Skipping."
        return
    fi

    # Check if source FASTA exists
    if [ ! -f "$FASTA_IN" ]; then
        echo "❌ ERROR: Custom pathogen FASTA file not found!"
        echo "   Please place your file at: $FASTA_IN"
        echo "   (You can change the expected filename at the top of this script)"
        exit 1
    fi
    
    echo "Building DIAMOND database from $FASTA_IN..."
    
    # Build the diamond database
    diamond makedb --in "$FASTA_IN" -d "$DB_OUT_BASE" --threads 4
    
    if [ -f "$DB_OUT_FINAL" ]; then
        echo "✅ Successfully created DIAMOND database: $DB_OUT_FINAL"
    else
        echo "❌ ERROR: DIAMOND database creation failed."
        exit 1
    fi
}

setup_swissprot_db() {
    echo -e "\n--- Setting up SwissProt + COG Combined Database ---"
    
    # Final combined database name
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