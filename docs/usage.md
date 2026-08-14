# Command-line usage

Run `metaquest COMMAND --help` for the authoritative option list.

## Global options

Global options appear before the command:

~~~bash
metaquest --version
metaquest --config metaquest.yaml analyze ...
metaquest --verbose analyze ...
metaquest --debug analyze ...
metaquest --quiet check
metaquest --no-color setup-db --list
~~~

| Option | Meaning |
|---|---|
| `--config FILE` | YAML configuration override |
| `--verbose` | additional progress context |
| `--debug` | commands and diagnostic output |
| `--quiet` | errors only |
| `--no-color` | disable ANSI colors and animation |

## Check

~~~bash
metaquest check --db-dir /data/metaquest-db
~~~

`check` validates required executables, Python packages, and database files.
The full check currently includes the provisional functional database.

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

## Analyze

### Taxonomy-only

~~~bash
metaquest analyze \
  --db-dir /data/metaquest-db \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --skip-annotation \
  --output results/sample
~~~

This runs validation, Kraken2, Bracken, descriptive reporting, and taxonomic
visualization.

### Full provisional path

~~~bash
metaquest analyze \
  --db-dir /data/metaquest-db \
  --single reads.fastq.gz \
  --output results/sample
~~~

This additionally runs MEGAHIT, Pyrodigal, and provisional DIAMOND annotation.
It requires a compatible functional database that the current database manager
does not yet install.

Useful options:

| Option | Meaning |
|---|---|
| `--output DIR` | result directory |
| `--db-dir DIR` | database root override |
| `--skip-validation` | bypass input validation |
| `--skip-annotation` | taxonomy-only execution |
| `--annotation-threads N` | DIAMOND thread allocation |

## Compare

~~~bash
metaquest compare \
  --inputs results/sample1 results/sample2 results/sample3 \
  --metadata metadata.tsv \
  --output comparison/
~~~

Minimum metadata:

~~~text
sample_id	group
sample1	control
sample2	control
sample3	treatment
~~~

The current comparison implementation is basic. Longitudinal models,
repeated-measures designs, ANCOM-BC2, and MaAsLin2 are planned but not yet part
of the maintained scientific contract.

## Setup databases

~~~bash
metaquest setup-db --db-dir /data/metaquest-db --list
metaquest setup-db --db-dir /data/metaquest-db --database taxonomy
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
├── functional_annotations.tsv
├── functional_annotations_filtered.tsv
├── 01_taxonomic_report.txt
├── 02_functional_report.txt
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
