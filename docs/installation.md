# Installation

This guide covers the maintained MetaQuest `2.0.0-alpha.1` runtime.

## Supported environment

- Linux or WSL2
- x86-64
- Python 3.10, 3.11, or 3.12
- conda or mamba
- at least 16 GB RAM recommended for the current Standard-8 profile
- at least 15 GB free during taxonomy database installation

MetaQuest has not yet been validated on native Windows or macOS.

## Install with conda

From the repository root:

~~~bash
conda env create -f environment/environment.yml
conda activate metaquest
python -m pip install -e .
~~~

For a regular non-editable installation:

~~~bash
python -m pip install .
~~~

Confirm the CLI:

~~~bash
metaquest --version
metaquest --help
~~~

## Runtime components

The minimal environment contains:

| Component | Role |
|---|---|
| Kraken2 | read-level taxonomic classification |
| Bracken | abundance re-estimation |
| MEGAHIT | metagenomic assembly |
| Pyrodigal | metagenomic gene prediction |
| eggNOG-mapper and DIAMOND | installed for the next functional milestone |
| BBMap | interleaved FASTQ splitting compatibility |

Pyrodigal is a Python dependency and does not require a separate Prodigal or
Prokka executable.

## Configure databases

Choose a database root outside Git:

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
metaquest databases --db-dir /data/metaquest-db --database taxonomy
metaquest databases --db-dir /data/metaquest-db --database functional
~~~

Then verify the installed tools and databases:

~~~bash
metaquest check --db-dir /data/metaquest-db
~~~

For repository-local convenience with physical storage on another disk:

~~~bash
mkdir -p /mnt/e/metaquest/databases
ln -s /mnt/e/metaquest/databases databases
metaquest databases --database taxonomy
~~~

See [Database management](databases.md) for manifests, checksums, and recovery.

## Configuration

Generate a YAML file:

~~~bash
metaquest init-config --output metaquest.yaml
~~~

Global options may appear before or after the command:

~~~bash
metaquest run --config metaquest.yaml \
  --single reads.fastq.gz \
  --skip-annotation \
  --output results/
~~~

## Development dependencies

Install the test extras:

~~~bash
python -m pip install -e '.[test]'
pytest -q
~~~

## Troubleshooting

### Database not found

Confirm the selected root:

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
ls -lh /data/metaquest-db/taxonomy/hash.k2d
~~~

Pass the same root to analysis or export it:

~~~bash
export METAQUEST_DB_DIR=/data/metaquest-db
~~~

### Interrupted database download

Run the same setup command again. Partial archives are retained under
`.downloads/` and resumed when the server supports HTTP range requests.

### Invalid existing database

MetaQuest will not silently overwrite an invalid directory. Inspect or move it,
or deliberately replace it:

~~~bash
metaquest databases \
  --db-dir /data/metaquest-db \
  --database taxonomy \
  --force
~~~

### Functional annotation is not present in analysis output

The legacy SwissProt/DIAMOND path has been removed. The installed eggNOG 5.0.2
database will be connected to the pipeline in the next functional milestone.
