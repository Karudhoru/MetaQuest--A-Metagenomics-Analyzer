"""Metagenomic gene prediction with Pyrodigal."""

import json
from pathlib import Path

from Bio import SeqIO

from ..exceptions import AnnotationError


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
