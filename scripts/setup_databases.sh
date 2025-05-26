#!/bin/bash
# Smart database setup that checks existing files
set -e

DB_DIR="$(dirname $0)/../databases"
mkdir -p "$DB_DIR"

# Check for existing Kraken2 database
if [ ! -f "$DB_DIR/hash.k2d" ] || [ ! -f "$DB_DIR/taxo.k2d" ]; then
    echo "Downloading MiniKraken database (8GB)"
    wget -nc -P "$DB_DIR" https://genome-idx.s3.amazonaws.com/kraken/minikraken_8GB_202003.tgz
    tar -xzf "$DB_DIR/minikraken_8GB_202003.tgz" -C "$DB_DIR"
else
    echo "Kraken2 database already exists - skipping download"
fi

# Pathogen database setup
PATHOGEN_DIR="$DB_DIR/pathogen_db"
mkdir -p "$PATHOGEN_DIR"

# CARD AMR Database
if [ ! -f "$PATHOGEN_DIR/nucleotide_fasta_protein_homolog_model.fasta" ]; then
    echo "Downloading CARD database..."
    wget -nc -P "$PATHOGEN_DIR" https://card.mcmaster.ca/latest/data 
    tar -xjf "$PATHOGEN_DIR/data" -C "$PATHOGEN_DIR"
else
    echo "CARD database exists - skipping download"
fi

# VFDB Virulence Factors
if [ ! -f "$PATHOGEN_DIR/VFDB_setB_pro.fas" ]; then
    echo "Downloading VFDB..."
    wget -nc -P "$PATHOGEN_DIR" http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas
else
    echo "VFDB exists - skipping download"
fi

# Build pathogen DIAMOND database if missing
if [ ! -f "$PATHOGEN_DIR/pathogen_db.dmnd" ]; then
    echo "Building pathogen DIAMOND database..."
    diamond makedb \
        --in "$PATHOGEN_DIR/nucleotide_fasta_protein_homolog_model.fasta" \
        --db "$PATHOGEN_DIR/pathogen_db"
else
    echo "Pathogen DIAMOND database exists - skipping creation"
fi

echo "Database setup complete!"
