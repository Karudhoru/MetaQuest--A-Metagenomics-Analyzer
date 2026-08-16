# Gene prediction and functional annotation

This document describes the maintained annotation boundary in
MetaQuest `2.0.0a1`.

## Current stage

~~~text
MEGAHIT contigs
  ↓
Pyrodigal metagenomic gene prediction
  ├── genes.faa
  ├── genes.fna
  ├── genes.gff3
  └── summary.json
  ↓
eggNOG-mapper 2.1.15 + eggNOG 5.0.2
  ├── functional_annotations.tsv
  ├── functional_category_summary.tsv
  ├── metaquest.emapper.annotations
  ├── summary.json
  └── completion.json
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

## eggNOG functional annotation

Every protein predicted by Pyrodigal is submitted to eggNOG-mapper using the
DIAMOND method and automatic taxonomic scope. MetaQuest left-joins mapper
results back to the complete gene catalog, so genes without an assignment
remain explicit `unannotated` rows rather than disappearing.

DIAMOND uses a `0.5` block size by default to keep the maintained workflow
practical on the supported 8-16 GB memory target. This value is recorded in
the functional summary and can be changed in YAML as
`annotation.diamond_block_size`.

| File | Contents |
|---|---|
| `functional_annotations.tsv` | one row per predicted gene, including unannotated genes |
| `functional_category_summary.tsv` | gene counts aggregated by COG, KO, EC, and GO term |
| `metaquest.emapper.annotations` | unmodified eggNOG-mapper annotation output |
| `eggnog_mapper.log` | mapper stdout and stderr for diagnostics |
| `summary.json` | tool/database versions, parameters, and annotation counts |
| `completion.json` | input checksum and parameters used for restart-safe reuse |

If the protein checksum, mapper version, database release, taxonomic scope, and
E-value match a completed run, MetaQuest reuses it. Partial or incompatible
outputs are replaced through eggNOG-mapper's explicit override mode.

These outputs describe gene presence, not abundance. Publication-ready
quantitative functional abundance still requires a validated read-mapping and
normalization stage.

## Interpretation limits

- A sequence-similarity hit is not experimental confirmation of function.
- Presence of a gene is not proof of expression or phenotype.
- Assembly can merge, fragment, or omit low-abundance sequences.
- Closely related proteins can have different substrate specificity.
- Functional abundance depends on read mapping and normalization choices.
