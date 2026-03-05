# MetaQuest Annotation Guide

This guide covers the functional annotation system used by MetaQuest, including database setup, annotation workflow, and interpreting results.

---

## Table of Contents

- [Overview](#overview)
- [Databases](#databases)
- [Annotation Workflow](#annotation-workflow)
- [Output Files](#output-files)
- [Interpreting Results](#interpreting-results)
- [Troubleshooting](#troubleshooting)

---

## Overview

MetaQuest uses a dual-database annotation strategy to maximise functional coverage:

| Database | Tool | Focus |
|----------|------|-------|
| **SwissProt** | DIAMOND | Manually reviewed protein annotations, pathogen markers |
| **COG** | DIAMOND | Clusters of Orthologous Genes — functional categories |

Gene prediction is performed by **Prokka** before annotation. DIAMOND then searches each predicted protein against both databases.

> **Note on annotation counts:** DIAMOND reports multiple database matches per query protein. One protein can match several SwissProt/COG entries. MetaQuest reports distinguish between *unique annotated proteins* (CDS-level) and *total annotation alignments* to avoid inflated counts.

---

## Databases

### SwissProt

- Source: [UniProt Knowledge Base](https://www.uniprot.org/uniprot/?query=reviewed:yes)
- Manually reviewed, high-quality protein annotations
- Covers virulence factors, AMR markers, housekeeping proteins
- Used for: functional annotation, pathogen marker detection

### COG (Clusters of Orthologous Genes)

- Source: [NCBI COG 2020](https://ftp.ncbi.nih.gov/pub/COG/COG2020/)
- Organized into 26 functional categories (e.g., J = Translation, K = Transcription, L = DNA replication)
- Used for: broad functional classification, metabolic pathway inference

### Combined DIAMOND Database

Both databases are merged into a single DIAMOND-formatted database (`SwissProt_COG_db.dmnd`) for a single-pass annotation run. This is built by `scripts/setup_databases.sh --swissprot`.

---

## Annotation Workflow

```
FASTQ/FASTA input
      │
      ▼
  Assembly (SPAdes/MEGAHIT)
      │
      ▼
  Gene Prediction (Prokka)
      │  → prokka_annotation/
      │     ├── *.faa   (protein sequences)
      │     ├── *.gff   (genome annotation)
      │     └── *.gbk   (GenBank format)
      ▼
  DIAMOND search vs SwissProt+COG
      │  → functional_annotations.tsv
      ▼
  Functional analysis & reporting
      │  → functional_report.txt
      └  → visualizations/
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-contig-length` | 500 bp | Minimum contig length for Prokka annotation |
| `--no-filter-contigs` | off | Annotate all contigs regardless of length |
| `--skip-annotation` | off | Skip Prokka+DIAMOND entirely (taxonomy only) |
| `--threads` | 4 | Threads for Prokka and DIAMOND |

---

## Output Files

| File | Description |
|------|-------------|
| `functional_annotations.tsv` | DIAMOND output: query, subject, identity, coverage, e-value, bitscore |
| `prokka_annotation/*.faa` | Predicted protein sequences |
| `prokka_annotation/*.gff` | Gene prediction coordinates |
| `functional_report.txt` | Summary with COG categories, annotation stats, methodology note |
| `cog_category_distribution.png` | Bar chart of COG functional categories |

### functional_annotations.tsv Columns

| Column | Description |
|--------|-------------|
| `query_id` | Prokka protein ID |
| `subject_id` | SwissProt/COG database hit |
| `pident` | Percentage identity |
| `length` | Alignment length |
| `qcovhsp` | Query coverage (%) |
| `evalue` | E-value |
| `bitscore` | Bit score |

---

## Interpreting Results

### Annotation Coverage

- **Annotated CDS**: Unique proteins with at least one database hit
- **Total annotations**: Sum of all DIAMOND alignments (≥ Annotated CDS, due to multi-hit behavior)
- A typical well-assembled bacterial genome achieves 85–95% annotation coverage against SwissProt+COG

### COG Functional Categories

| Letter | Category |
|--------|----------|
| J | Translation, ribosomal structure |
| K | Transcription |
| L | DNA replication, recombination, repair |
| D | Cell cycle, division |
| T | Signal transduction |
| M | Cell wall, membrane, envelope |
| N | Cell motility |
| O | Post-translational modification |
| X | Mobilome (transposons, IS elements) |
| C | Energy metabolism |
| G | Carbohydrate metabolism |
| E | Amino acid metabolism |
| F | Nucleotide metabolism |
| H | Coenzyme metabolism |
| I | Lipid metabolism |
| P | Inorganic ion transport |
| Q | Secondary metabolite biosynthesis |
| U | Intracellular trafficking |
| V | **Defense mechanisms (AMR)** |
| R | General function prediction |
| S | Function unknown |

### Mobile Genetic Elements

MetaQuest tracks COG category **X** (Mobilome) separately, reporting:
- IS family classification
- Transposase count and diversity
- Integron and phage element signals

---

## Troubleshooting

**No annotations produced:**
- Verify database exists: `ls databases/SwissProt_COG_db.dmnd`
- Check Prokka ran successfully: `ls results/prokka_annotation/*.faa`
- Run with `--debug` to see DIAMOND command and output

**Very low coverage (<50%):**
- Assembly quality may be poor — check N50 in the taxonomic report
- Try lowering the DIAMOND e-value threshold in `config.py`

**Annotation takes too long:**
- Increase `--threads`
- Use `--min-contig-length 1000` to reduce the number of short contigs annotated

**Rebuilding the database:**
```bash
rm databases/SwissProt_COG_db.dmnd
./scripts/setup_databases.sh --swissprot
```
