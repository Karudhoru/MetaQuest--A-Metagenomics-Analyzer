# Command-line usage

Run `metaquest COMMAND --help` for the authoritative option list.

## Global options

Global options may appear before or after the command:

~~~bash
metaquest --version
metaquest run --config metaquest.yaml ...
metaquest run --verbose ...
metaquest run --debug ...
metaquest --quiet check
metaquest databases --no-color --list
~~~

| Option | Meaning |
|---|---|
| `--config FILE` | YAML configuration override |
| `--verbose` | additional progress context |
| `--debug` | commands and diagnostic output |
| `--quiet` | errors only |
| `--no-color` | disable ANSI colors and animation |
| `--low-memory` | use Kraken2 database memory mapping |

## Check

~~~bash
metaquest check --db-dir /data/metaquest-db
~~~

`check` validates required executables, Python packages, and database files.
The full check validates the installed eggNOG core database.

## Validate FASTQ

~~~bash
metaquest validate --single reads.fastq.gz

metaquest validate \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz

metaquest validate --interleaved sample.fastq.gz
~~~

Optional thresholds:

~~~text
--min-quality SCORE
--min-sequences NUM
--overrep-threshold FRACTION
~~~

Validation does not currently trim reads. fastp integration is planned for the
Snakemake workflow.

## Run

### Taxonomy-only

~~~bash
metaquest run \
  --db-dir /data/metaquest-db \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --taxonomy-only \
  --output results/sample
~~~

This runs validation, Kraken2, Bracken, and descriptive reporting.

### Assembly, gene prediction, and functional annotation

~~~bash
metaquest run \
  --db-dir /data/metaquest-db \
  --single reads.fastq.gz \
  --output results/sample
~~~

This additionally runs MEGAHIT, Pyrodigal, and eggNOG-mapper. Use
`--skip-functional` when only the predicted gene catalog is required.

Useful options:

| Option | Meaning |
|---|---|
| `--output DIR` | result directory |
| `--db-dir DIR` | database root override |
| `--skip-validation` | bypass input validation |
| `--taxonomy-only` | taxonomy-only execution |
| `--skip-annotation` | deprecated alias for `--taxonomy-only` |
| `--skip-functional` | stop after Pyrodigal gene prediction |
| `--annotation-threads N` | eggNOG-mapper worker threads |

## Compare

The `compare` command name is retained for compatibility, but the experimental
legacy implementation has been removed. It exits with an explanatory error
until a validated comparative workflow is built.

## Setup databases

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
metaquest databases --db-dir /data/metaquest-db --database taxonomy
~~~

Use `--force` only when deliberately replacing an existing invalid or older
installation.

See [Database management](databases.md).

## Initialize configuration

~~~bash
metaquest init-config --output metaquest.yaml
~~~

Database root resolution:

~~~text
--db-dir
  → METAQUEST_DB_DIR
    → databases.base_dir in YAML
      → ./databases
~~~

## Current output layout

~~~text
results/sample/
├── kraken_classified.txt
├── kraken_report.txt
├── bracken_report.tsv
├── megahit_assembly/
├── gene_prediction/
│   ├── genes.faa
│   ├── genes.fna
│   ├── genes.gff3
│   └── summary.json
├── 01_taxonomic_report.txt
├── analysis_summary.json
├── analysis_metadata.json
└── metaquest.log
~~~

Taxonomy-only runs omit assembly and annotation outputs.

## Exit behavior

- `0`: command completed successfully
- `1`: validation, dependency, database, or pipeline failure
- `130`: interrupted with Ctrl+C

Errors include an actionable message. Detailed external-tool output is retained
in the run log when logging is enabled.
