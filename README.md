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
- descriptive text and JSON reporting
- versioned taxonomy and eggNOG database installation

MetaQuest does not currently perform read trimming, host-read removal, validated
AMR or virulence analysis, clinical diagnosis, or treatment recommendations.

## Pipeline

~~~text
FASTQ
  ├── validation
  ├── Kraken2 → Bracken
  ├── MEGAHIT
  ├── Pyrodigal
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
metaquest databases --list
metaquest databases --database taxonomy
metaquest check
~~~

By default, MetaQuest installs reference data in `./databases`. Use
`--db-dir` or `METAQUEST_DB_DIR` only when you intentionally want another
location.

Run taxonomic profiling:

~~~bash
metaquest run \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --skip-annotation \
  --output results/sample
~~~

The installed eggNOG database is not yet connected to the analysis pipeline.
Functional annotation is the next core implementation milestone.

## Commands

~~~text
metaquest check       # verify tools and databases
metaquest validate    # validate FASTQ input
metaquest run         # run the analysis pipeline
metaquest databases   # inspect or install reference data
metaquest init-config # create a configuration file
~~~

The previous `analyze` and `setup-db` names remain available as aliases for
backward compatibility.

Global output controls may appear before or after the command:

~~~bash
metaquest run --verbose --single reads.fastq.gz --output results/
metaquest databases --no-color --list
metaquest --quiet check
~~~

For systems that cannot load the Kraken2 database fully into RAM, enable
Kraken2 memory mapping with `--low-memory`:

~~~bash
metaquest run --low-memory --single reads.fastq.gz --output results/
~~~

Currently this flag only adds Kraken2's `--memory-mapping` option. It does not
change resource settings for Bracken, MEGAHIT, or gene prediction.

## Database storage

Reference databases are not stored in Git. From the repository root, install
them directly into the default `databases/` directory:

~~~bash
metaquest databases --database taxonomy
metaquest databases --database functional
metaquest check
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
