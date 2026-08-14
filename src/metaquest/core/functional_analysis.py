"""Gene prediction and provisional DIAMOND functional annotation."""

import json
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from ..exceptions import AnnotationError
from ..io.output_formatter import get_formatter


def run_gene_prediction(
    fasta_path: Path,
    output_dir: Path,
    *,
    min_contig_length: int = 200,
) -> tuple[Path, int]:
    """Predict genes from metagenomic contigs with bounded-memory Pyrodigal."""
    try:
        import pyrodigal
    except ImportError as exc:
        raise AnnotationError(
            "Pyrodigal is required for gene prediction; install pyrodigal>=3.7.1"
        ) from exc

    prediction_dir = Path(output_dir) / "gene_prediction"
    prediction_dir.mkdir(parents=True, exist_ok=True)
    proteins_path = prediction_dir / "genes.faa"
    genes_path = prediction_dir / "genes.fna"
    gff_path = prediction_dir / "genes.gff3"

    finder = pyrodigal.GeneFinder(meta=True)
    contigs_seen = 0
    contigs_processed = 0
    genes_predicted = 0

    try:
        with proteins_path.open("w", encoding="utf-8") as proteins_out, \
             genes_path.open("w", encoding="utf-8") as genes_out, \
             gff_path.open("w", encoding="utf-8") as gff_out:
            gff_out.write("##gff-version 3\n")
            for record in SeqIO.parse(Path(fasta_path), "fasta"):
                contigs_seen += 1
                if len(record.seq) < min_contig_length:
                    continue
                predictions = finder.find_genes(bytes(record.seq))
                contigs_processed += 1
                genes_predicted += len(predictions)
                predictions.write_translations(
                    proteins_out,
                    sequence_id=record.id,
                    include_stop=False,
                    full_id=True,
                )
                predictions.write_genes(
                    genes_out,
                    sequence_id=record.id,
                    full_id=True,
                )
                predictions.write_gff(
                    gff_out,
                    sequence_id=record.id,
                    header=False,
                    include_translation_table=True,
                    full_id=True,
                )
    except Exception as exc:
        raise AnnotationError(f"Pyrodigal gene prediction failed: {exc}") from exc

    summary = {
        "tool": "Pyrodigal",
        "tool_version": getattr(pyrodigal, "__version__", "unknown"),
        "mode": "metagenomic",
        "minimum_contig_length": min_contig_length,
        "contigs_seen": contigs_seen,
        "contigs_processed": contigs_processed,
        "genes_predicted": genes_predicted,
    }
    (prediction_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return prediction_dir, genes_predicted


def run_functional_annotation(
    gene_prediction_dir: Path,
    output_dir: Path,
    *,
    threads: int = 8,
    evalue: float = 1e-5,
    sensitivity: str = "sensitive",
    db_path: Path | None = None,
) -> Optional[Path]:
    """
    Run DIAMOND functional annotation against SwissProt/COG.

    Args:
        gene_prediction_dir: Pyrodigal output directory containing genes.faa.
        output_dir: Output directory for annotation results.
        threads: Number of threads.
        evalue: E-value threshold.
        sensitivity: DIAMOND sensitivity mode.
        db_path: Path to DIAMOND database. Defaults to config.

    Returns:
        Path to annotation TSV file, or None on failure.
    """
    formatter = get_formatter()

    if db_path is None:
        from ..settings import get_config
        db_path = get_config().databases.swissprot_cog

    protein_fasta = Path(gene_prediction_dir) / "genes.faa"
    diamond_output = Path(output_dir) / "functional_annotations.tsv"

    if not protein_fasta.exists() or protein_fasta.stat().st_size == 0:
        formatter.warning("Protein FASTA missing or empty")
        return None

    if not db_path.exists():
        raise AnnotationError(f"SwissProt database not found: {db_path}")

    protein_count = sum(1 for _ in SeqIO.parse(protein_fasta, "fasta"))
    formatter.debug(f"Annotating {protein_count:,} proteins")

    outfmt_cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
    cmd = [
        "diamond", "blastp",
        "--db", str(db_path),
        "--query", str(protein_fasta),
        "--out", str(diamond_output),
        "--outfmt", "6", *outfmt_cols.split(),
        "--top", "1",
        "--evalue", str(evalue),
        "--threads", str(threads),
        f"--{sensitivity}",
        "--block-size", "4.0",
        "--index-chunks", "1",
        "--log",
    ]

    return_code, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="DIAMOND Annotation (SwissProt)",
        capture_output=True,
        show_command=False,
    )

    if return_code != 0:
        formatter.warning(f"DIAMOND failed (exit {return_code})")
        return None

    if not diamond_output.exists() or diamond_output.stat().st_size == 0:
        formatter.warning("DIAMOND produced no output")
        return None

    # Post-filter: remove low-quality hits
    filtered_output = Path(output_dir) / "functional_annotations_filtered.tsv"
    unique_proteins = set()
    removed = 0

    with open(diamond_output) as fin, open(filtered_output, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 15:
                fout.write(line)
                unique_proteins.add(parts[0])
                continue

            pident = float(parts[2])
            length = int(parts[3])
            qlen = int(parts[12])

            # Filter: min 40% identity AND min 50% query coverage
            query_coverage = length / qlen if qlen > 0 else 0
            if pident < 40.0 or query_coverage < 0.5:
                removed += 1
                continue

            fout.write(line)
            unique_proteins.add(parts[0])

    if removed > 0:
        formatter.debug(f"Post-filter removed {removed} low-quality annotations (<40% identity or <50% query coverage)")

    # Replace original with filtered
    filtered_output.replace(diamond_output)

    if not unique_proteins:
        formatter.warning("DIAMOND output contains no valid annotations after filtering")
        return None

    annotation_pct = (len(unique_proteins) / protein_count * 100) if protein_count > 0 else 0
    formatter.debug(f"Annotated {len(unique_proteins)}/{protein_count} proteins ({annotation_pct:.1f}%)")

    return diamond_output
