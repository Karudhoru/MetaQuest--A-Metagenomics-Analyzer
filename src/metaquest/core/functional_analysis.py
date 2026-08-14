"""
Functional Analysis Module — Prokka gene prediction + DIAMOND annotation.
"""

import subprocess
import os
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from ..exceptions import AnnotationError
from ..io.output_formatter import get_formatter


def parse_prokka_sample_txt(sample_file: Path) -> dict:
    """Parse Prokka's sample.txt file for gene counts."""
    stats = {
        "organism": "",
        "contigs": 0,
        "bases": 0,
        "CDS": 0,
        "cds": 0,
        "gene": 0,
        "rRNA": 0,
        "tRNA": 0,
        "repeat_region": 0,
    }

    if not Path(sample_file).exists():
        return stats

    with open(sample_file) as f:
        for line in f:
            line = line.strip()
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                key_lower = key.lower()

                if key_lower == "organism":
                    stats["organism"] = value
                elif key_lower in ("cds", "gene", "rrna", "trna", "contigs", "bases", "repeat_region"):
                    try:
                        if key_lower == "cds":
                            stats["CDS"] = int(value)
                            stats["cds"] = int(value)
                        elif key_lower == "rrna":
                            stats["rRNA"] = int(value)
                        elif key_lower == "trna":
                            stats["tRNA"] = int(value)
                        elif key_lower == "gene":
                            stats["gene"] = int(value)
                        else:
                            stats[key_lower] = int(value)
                    except ValueError:
                        pass

    return stats


def run_prokka(
    fasta_path: Path,
    output_dir: Path,
    *,
    kill_tbl2asn: bool = False,
    tbl2asn_timeout: int = 30,
    threads: int = 8,
    mode: str = "metagenome",
) -> Path:
    """
    Run Prokka gene prediction.

    Args:
        fasta_path: Input contigs FASTA.
        output_dir: Output directory.
        kill_tbl2asn: Deprecated compatibility argument; ignored.
        tbl2asn_timeout: Deprecated compatibility argument; ignored.
        threads: CPU threads.
        mode: Prokka mode (metagenome, careful, fast).

    Returns:
        Path to Prokka output directory.

    Raises:
        AnnotationError: If Prokka fails to produce essential output.
    """
    formatter = get_formatter()
    prokka_dir = output_dir / "prokka_annotation"

    cmd = [
        "prokka",
        "--outdir", str(prokka_dir),
        "--prefix", "sample",
        "--cpus", str(threads),
        "--centre", "X",
        "--compliant",
        str(fasta_path),
    ]

    if mode == "metagenome":
        cmd.insert(-1, "--metagenome")
        cmd.insert(-1, "--fast")
    elif mode == "careful":
        pass  # Prokka default is careful
    elif mode == "fast":
        cmd.insert(-1, "--fast")

    formatter.debug(f"Prokka command: {' '.join(cmd)}")

    return_code, _, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Prokka gene prediction",
        capture_output=True,
        show_command=False,
    )
    if return_code != 0:
        raise AnnotationError(f"Prokka failed (exit {return_code}): {stderr}")

    # Validate essential output
    essential = [prokka_dir / "sample.faa", prokka_dir / "sample.ffn", prokka_dir / "sample.gff"]
    missing = [f.name for f in essential if not f.exists()]
    if missing:
        raise AnnotationError(f"Prokka output is incomplete: {', '.join(missing)}")

    return prokka_dir


def run_functional_annotation(
    prokka_dir: Path,
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
        prokka_dir: Prokka output directory containing sample.faa.
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

    protein_fasta = Path(prokka_dir) / "sample.faa"
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
