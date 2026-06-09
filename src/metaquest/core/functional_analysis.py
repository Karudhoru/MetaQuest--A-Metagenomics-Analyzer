"""
Functional Analysis Module — Prokka gene prediction + DIAMOND annotation.
"""

import subprocess
import os
import time
import threading
import queue
from pathlib import Path
from typing import Optional

import psutil
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
    kill_tbl2asn: bool = True,
    tbl2asn_timeout: int = 30,
    threads: int = 8,
    mode: str = "metagenome",
) -> Path:
    """
    Run Prokka gene prediction.

    Args:
        fasta_path: Input contigs FASTA.
        output_dir: Output directory.
        kill_tbl2asn: Monitor and kill long-running tbl2asn processes.
        tbl2asn_timeout: Seconds before killing tbl2asn.
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
        "--force",
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

    prokka_process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Drain pipes in background threads to prevent deadlock
    stdout_queue: queue.Queue = queue.Queue()
    stderr_queue: queue.Queue = queue.Queue()

    def read_pipe(pipe, q):
        try:
            for line in iter(pipe.readline, ""):
                if line:
                    q.put(line)
            pipe.close()
        except Exception:
            pass

    stdout_thread = threading.Thread(target=read_pipe, args=(prokka_process.stdout, stdout_queue), daemon=True)
    stderr_thread = threading.Thread(target=read_pipe, args=(prokka_process.stderr, stderr_queue), daemon=True)
    stdout_thread.start()
    stderr_thread.start()

    # Monitor and kill tbl2asn if needed
    killed_count = 0
    if kill_tbl2asn:
        tbl2asn_first_seen = None
        while prokka_process.poll() is None:
            try:
                for proc in psutil.process_iter(["pid", "name"]):
                    try:
                        if "tbl2asn" in proc.info.get("name", "").lower():
                            if tbl2asn_first_seen is None:
                                tbl2asn_first_seen = time.time()
                            if time.time() - tbl2asn_first_seen > tbl2asn_timeout:
                                proc.kill()
                                proc.wait(timeout=2)
                                killed_count += 1
                                tbl2asn_first_seen = None
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue
            except Exception:
                pass
            time.sleep(0.5)
    else:
        prokka_process.wait()

    stdout_thread.join(timeout=5)
    stderr_thread.join(timeout=5)

    if killed_count > 0:
        formatter.warning(f"Killed {killed_count} tbl2asn process(es)")

    # Cleanup any remaining tbl2asn
    if kill_tbl2asn:
        try:
            subprocess.run(["pkill", "-9", "tbl2asn"], stderr=subprocess.DEVNULL, timeout=2)
        except Exception:
            pass

    # Validate essential output
    essential = [prokka_dir / "sample.faa", prokka_dir / "sample.ffn", prokka_dir / "sample.gff"]
    missing = [f.name for f in essential if not f.exists()]
    if missing:
        formatter.warning(f"Missing Prokka files: {', '.join(missing)}")

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

    outfmt_cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
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

    # Validate
    unique_proteins = set()
    with open(diamond_output) as f:
        for line in f:
            if line.strip():
                unique_proteins.add(line.split("\t")[0])

    if not unique_proteins:
        formatter.warning("DIAMOND output contains no valid annotations")
        return None

    annotation_pct = (len(unique_proteins) / protein_count * 100) if protein_count > 0 else 0
    formatter.debug(f"Annotated {len(unique_proteins)}/{protein_count} proteins ({annotation_pct:.1f}%)")

    return diamond_output
