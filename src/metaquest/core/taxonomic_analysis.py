"""
Taxonomic Analysis Module — Kraken2/Bracken
"""

import gzip
import re
import time
from collections import Counter
from pathlib import Path

from ..exceptions import ClassificationError
from ..io.output_formatter import get_formatter


def infer_read_length(input_files, *, sample_reads: int = 10_000) -> int:
    """Return the modal read length from a bounded FASTQ sample."""
    lengths: Counter[int] = Counter()
    remaining = sample_reads
    for input_file in input_files:
        opener = gzip.open if str(input_file).endswith(".gz") else open
        with opener(input_file, "rt", encoding="utf-8") as handle:
            while remaining:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip("\r\n")
                plus = handle.readline()
                quality = handle.readline().rstrip("\r\n")
                if (
                    not header.startswith("@")
                    or not plus.startswith("+")
                    or len(sequence) != len(quality)
                ):
                    raise ClassificationError(
                        f"Invalid FASTQ record while inferring read length: {input_file}"
                    )
                lengths[len(sequence)] += 1
                remaining -= 1
            if remaining == 0:
                break
    if not lengths:
        raise ClassificationError("Cannot infer read length from empty FASTQ input")
    return lengths.most_common(1)[0][0]


def select_bracken_read_length(db_path: Path, observed: int) -> int:
    """Select the closest installed Bracken k-mer distribution model."""
    available = []
    for path in Path(db_path).glob("database*mers.kmer_distrib"):
        match = re.fullmatch(r"database(\d+)mers\.kmer_distrib", path.name)
        if match:
            available.append(int(match.group(1)))
    if not available:
        raise ClassificationError(f"No Bracken read-length models found in {db_path}")
    return min(available, key=lambda value: (abs(value - observed), value))


def parse_kraken_read_counts(report_path: Path) -> tuple[int, int, int]:
    """Return total, classified, and unclassified reads from a Kraken report."""
    unclassified = 0
    root_reads = 0
    with Path(report_path).open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            rank = fields[3].strip()
            if rank == "U":
                unclassified = int(fields[1])
            elif rank == "R":
                root_reads = int(fields[1])
    total = root_reads + unclassified
    return total, root_reads, unclassified


def run_kraken(
    input_files,
    output_dir,
    *,
    threads: int = 8,
    min_hit_groups: int = 2,
    db_path: Path | None = None,
    memory_mapping: bool = False,
):
    """
    Run Kraken2 classification for FASTQ files.

    Args:
        input_files: FASTQ file path(s) — single file or list of paired files.
        output_dir: Output directory for results.
        threads: Number of threads.
        min_hit_groups: Minimum hit groups for classification.
        db_path: Kraken2 database path. Defaults to config value.
        memory_mapping: Use Kraken2 memory mapping instead of loading the database into RAM.

    Returns:
        Path to Kraken2 report file.
    """
    formatter = get_formatter()

    if db_path is None:
        from ..settings import get_config
        db_path = get_config().databases.kraken_db

    report = Path(output_dir) / "kraken_report.txt"
    classified = Path(output_dir) / "kraken_classified.txt"

    cmd = [
        "kraken2",
        "--db", str(db_path),
        "--threads", str(threads),
        "--minimum-hit-groups", str(min_hit_groups),
        "--report", str(report),
        "--output", str(classified),
    ]

    if memory_mapping:
        cmd.append("--memory-mapping")

    if isinstance(input_files, (list, tuple)) and len(input_files) == 2:
        cmd.extend(["--paired", str(input_files[0]), str(input_files[1])])
    else:
        input_file = input_files[0] if isinstance(input_files, (list, tuple)) else input_files
        cmd.append(str(input_file))

    formatter.debug(f"Kraken2 command: {' '.join(cmd)}")

    returncode, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Kraken2 classification",
        capture_output=True,
        show_command=False,
    )

    if returncode != 0:
        raise ClassificationError(f"Kraken2 failed (exit {returncode}): {stderr}")

    if not report.exists() or report.stat().st_size == 0:
        raise ClassificationError(f"Kraken2 report missing or empty: {report}")

    return report


def run_bracken(
    report_path,
    output_dir,
    *,
    read_length: int = 150,
    level: str = "S",
    threshold: int = 10,
    db_path: Path | None = None,
):
    """
    Estimate species abundances with Bracken.

    Args:
        report_path: Path to Kraken2 report file.
        output_dir: Output directory.
        read_length: Sequencing read length.
        level: Taxonomic level (S=Species, G=Genus, F=Family).
        threshold: Minimum read threshold.
        db_path: Kraken database path.

    Returns:
        Path to the Bracken output file.
    """
    formatter = get_formatter()

    if db_path is None:
        from ..settings import get_config
        db_path = get_config().databases.kraken_db

    bracken_out = Path(output_dir) / "bracken_report.tsv"
    report_file = Path(report_path)

    # Wait for Kraken report
    for _ in range(5):
        if report_file.exists() and report_file.stat().st_size > 0:
            break
        time.sleep(1)
    else:
        raise ClassificationError(f"Kraken2 report not found: {report_file}")

    cmd = [
        "bracken",
        "-d", str(db_path),
        "-i", str(report_path),
        "-o", str(bracken_out),
        "-w", str(Path(output_dir) / "bracken_report.txt"),
        "-r", str(read_length),
        "-l", level,
        "-t", str(threshold),
    ]

    formatter.debug(f"Bracken: level={level}, read_len={read_length}, threshold={threshold}")

    returncode, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Bracken abundance estimation",
        capture_output=True,
        show_command=False,
    )

    if returncode != 0:
        raise ClassificationError(
            f"Bracken failed (exit {returncode}): {stderr}"
        )

    if not bracken_out.exists() or bracken_out.stat().st_size == 0:
        raise ClassificationError(f"Bracken report missing or empty: {bracken_out}")

    return bracken_out
