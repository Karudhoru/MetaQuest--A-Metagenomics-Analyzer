# MetaQuest

MetaQuest is a research-use command-line pipeline for short-read metagenomic
FASTQ analysis.

The project is undergoing a production-readiness refactor. Version
`2.0.0-alpha.1`
deliberately exposes a smaller stable surface while the final QC, functional,
pathogen-associated, database, and comparative workflows are evaluated.

## Current stable surface

- FASTQ input and quality validation
- Kraken2 taxonomic classification
- Bracken abundance estimation
- MEGAHIT metagenomic assembly
- descriptive taxonomic and assembly outputs
- provisional Prokka + DIAMOND functional annotation
- basic taxonomic group comparison
- text, JSON, TSV, and taxonomic visualization outputs

MetaQuest is for scientific research. It does not provide clinical diagnosis,
patient risk assessment, or treatment recommendations.

## Experimental features

The repository contains research prototypes for custom pathogen markers,
protein machine learning, HMM profiles, ESM embeddings, pathogenicity-island
detection, and Bayesian risk integration. These features are not part of the
default pipeline and are disabled until independently validated.

Their presence in the source tree must not be interpreted as evidence of
accuracy, clinical utility, or publication readiness.

## Commands

The existing command structure is retained during the refactor:

```text
metaquest check
metaquest validate
metaquest analyze
metaquest compare
metaquest setup-db
metaquest init-config
```

### Validate FASTQ input

```bash
metaquest validate --single reads.fastq.gz
metaquest validate --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

Validation currently assesses the input; it does not trim or clean reads.
The production QC tool will be selected in the next workflow-design phase.

### Taxonomic-only analysis

```bash
metaquest analyze \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --skip-annotation \
  --output results/
```

This executes Kraken2, Bracken, descriptive reporting, and taxonomic
visualization. Assembly and annotation are skipped.

### Current full analysis

```bash
metaquest analyze \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --output results/
```

The current full path adds MEGAHIT and provisional Prokka/DIAMOND annotation.
Prokka and the functional database strategy are temporary pending the tool and
database review. No custom pathogen or risk model runs in this workflow.

### Basic comparison

```bash
metaquest compare \
  --inputs sample1/ sample2/ sample3/ sample4/ \
  --metadata metadata.tsv \
  --output comparison/
```

The current comparison accepts a `sample_id` and `group` design. Longitudinal,
paired, repeated-measures, functional, AMR, and virulence comparisons are
planned but are not currently claimed as implemented.

## Output

Stable analysis outputs include:

```text
results/
├── kraken_report.txt
├── bracken_report.tsv
├── 01_taxonomic_report.txt
├── analysis_summary.json
├── analysis_metadata.json
├── megahit_assembly/              # full path only
├── prokka_annotation/             # provisional full path only
├── functional_annotations.tsv     # provisional full path only
└── 02_functional_report.txt       # provisional full path only
```

`analysis_summary.json` explicitly records that experimental features were not
executed.

## Configuration

Generate an editable configuration file with:

```bash
metaquest init-config --output metaquest.yaml
```

Use it with:

```bash
metaquest --config metaquest.yaml analyze --single reads.fastq.gz
```

Database paths can also be rooted with `METAQUEST_DB_DIR`.

## Development status

The immediate priorities are:

1. stabilize and test the CLI and current stage contracts;
2. select the production QC, functional, and pathogen-associated tools;
3. design metadata-aware condition and time-series comparison;
4. establish reproducible benchmark datasets and database manifests;
5. validate every scientific claim before publication.

The detailed stabilization plan is available in
[`docs/METAQUEST_ROADMAP.md`](docs/METAQUEST_ROADMAP.md).

## Installation

MetaQuest currently targets Linux and conda-based bioinformatics environments.
The existing environment file is retained during tool selection, but it is not
yet the final minimal production environment.

See [`docs/installation.md`](docs/installation.md) for the legacy installation
notes. Those notes are being revised as part of the refactor.

## License and citation

The repository declares the MIT license in its package metadata. Formal
citation information will be added after the workflow and benchmark are
published.
