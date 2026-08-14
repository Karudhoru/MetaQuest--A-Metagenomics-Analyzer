# MetaQuest

MetaQuest is a research-use command-line pipeline for short-read metagenomic
FASTQ analysis. Version `2.0.0-alpha.1` is a stabilization release: the
maintained runtime is intentionally smaller than earlier prototypes and does
not make clinical or pathogen-risk claims.

## Current capabilities

- FASTQ format and basic quality validation
- Kraken2 taxonomic classification
- Bracken abundance re-estimation
- optional taxonomy-only execution
- MEGAHIT metagenomic assembly
- Pyrodigal metagenomic gene prediction
- provisional DIAMOND functional similarity search
- basic metadata-driven group comparison
- descriptive text, JSON, TSV, HTML, and plot outputs
- versioned, checksummed taxonomy database installation

MetaQuest does not currently perform read trimming, host-read removal, validated
AMR or virulence analysis, clinical diagnosis, or treatment recommendations.

## Pipeline

~~~text
FASTQ
  ├── validation
  ├── Kraken2 → Bracken
  ├── MEGAHIT
  ├── Pyrodigal
  ├── provisional DIAMOND annotation
  └── descriptive reporting
~~~

The maintained pipeline currently runs directly in Python. Migration to
Snakemake, fastp QC, and optional Bowtie2 host filtering is planned but is not
claimed as implemented.

## Quick start

Create the environment and install MetaQuest:

~~~bash
conda env create -f environment/environment.yml
conda activate metaquest
python -m pip install -e .
~~~

Inspect and install the taxonomy database:

~~~bash
metaquest setup-db --db-dir /data/metaquest-db --list
metaquest setup-db --db-dir /data/metaquest-db --database taxonomy
metaquest check --db-dir /data/metaquest-db
~~~

Run taxonomic profiling:

~~~bash
metaquest analyze \
  --db-dir /data/metaquest-db \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --skip-annotation \
  --output results/sample
~~~

The complete non-taxonomy path remains provisional until the functional
database and abundance contract are finalized.

## Commands

~~~text
metaquest check
metaquest validate
metaquest analyze
metaquest compare
metaquest setup-db
metaquest init-config
~~~

Global output controls must appear before the command:

~~~bash
metaquest --verbose analyze --single reads.fastq.gz --output results/
metaquest --no-color setup-db --list
metaquest --quiet check
~~~

## Database storage

Reference databases are not stored in Git. A repository-local `databases`
path may be linked to a larger disk:

~~~bash
mkdir -p /mnt/e/metaquest/databases
ln -s /mnt/e/metaquest/databases databases
metaquest setup-db --database taxonomy
~~~

The database path resolution order is:

1. `--db-dir`
2. `METAQUEST_DB_DIR`
3. `databases.base_dir` in YAML
4. `./databases`

## Documentation

- [Installation](docs/installation.md)
- [Database management](docs/databases.md)
- [Usage](docs/usage.md)
- [Gene prediction and functional annotation](docs/annotation.md)
- [Changelog](CHANGELOG.md)

## Development status

MetaQuest is alpha software intended for reproducible method development. Before
a stable release, the workflow requires end-to-end Snakemake execution,
functional and pathogen-associated method validation, resource benchmarks,
comparison-method validation, continuous integration, and publication datasets.

## License

MetaQuest is distributed under the MIT license. Some runtime dependencies use
other licenses; users and redistributors must comply with each dependency's
terms. Pyrodigal is distributed under GPL-3.0-or-later.
