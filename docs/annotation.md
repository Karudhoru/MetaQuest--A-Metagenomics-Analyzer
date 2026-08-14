# Gene prediction and functional annotation

This document describes the maintained annotation boundary in
MetaQuest `2.0.0-alpha.1`.

## Current stage

~~~text
MEGAHIT contigs
  ↓
Pyrodigal metagenomic gene prediction
  ├── genes.faa
  ├── genes.fna
  ├── genes.gff3
  └── summary.json
~~~

## Pyrodigal

MetaQuest initializes `pyrodigal.GeneFinder(meta=True)` and processes one
contig at a time. This bounds memory use and applies Pyrodigal's metagenomic
models independently to assembled contigs.

Contigs shorter than `annotation.min_contig_length` are skipped. The default is
200 bp.

### Outputs

| File | Contents |
|---|---|
| `gene_prediction/genes.faa` | predicted protein sequences without terminal stop characters |
| `gene_prediction/genes.fna` | predicted coding nucleotide sequences |
| `gene_prediction/genes.gff3` | CDS coordinates and translation-table metadata |
| `gene_prediction/summary.json` | tool version, mode, contig counts, and gene count |

Sequence identifiers derive from input contig identifiers. Input contig IDs
therefore need to be unique.

## Planned production functional contract

The intended implementation is:

~~~text
cleaned reads ───────────────────────────────┐
                                            ↓
contigs → Pyrodigal genes → eggNOG annotation
                         ↑                  ↓
                         └── read mapping ──┘
                                            ↓
                           gene and function abundance
~~~

Publication-ready functional analysis requires both annotation and quantitative
read support. Planned output levels include raw gene counts and aggregated COG,
KO, EC, GO, and pathway tables. Database version, mapper parameters,
normalization, and aggregation rules must be recorded in provenance.

## Interpretation limits

- A sequence-similarity hit is not experimental confirmation of function.
- Presence of a gene is not proof of expression or phenotype.
- Assembly can merge, fragment, or omit low-abundance sequences.
- Closely related proteins can have different substrate specificity.
- Functional abundance depends on read mapping and normalization choices.
