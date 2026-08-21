"""fastp read preprocessing and machine-readable QC parsing."""

from __future__ import annotations

import json
from pathlib import Path

from metaquest.exceptions import MetaQuestError
from metaquest.io.output_formatter import get_formatter


def _summary_metrics(payload: dict) -> dict:
    """Extract stable, report-oriented values from fastp JSON."""
    summary = payload.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = payload.get("filtering_result", {})

    def metrics(values):
        total_bases = int(values.get("total_bases", 0) or 0)
        read1_length = float(values.get("read1_mean_length", 0.0) or 0.0)
        read2_value = values.get("read2_mean_length")
        read2_length = float(read2_value or 0.0) if read2_value is not None else None
        return {
            "reads": int(values.get("total_reads", 0) or 0),
            "bases": total_bases,
            "q20_bases": int(values.get("q20_bases", 0) or 0),
            "q30_bases": int(values.get("q30_bases", 0) or 0),
            "q20_rate": float(values.get("q20_rate", 0.0) or 0.0),
            "q30_rate": float(values.get("q30_rate", 0.0) or 0.0),
            "gc_content": float(values.get("gc_content", 0.0) or 0.0),
            "read1_mean_length": read1_length,
            "read2_mean_length": read2_length,
            "mean_length": (
                (read1_length + read2_length) / 2
                if read2_length is not None
                else read1_length
            ),
        }

    before_metrics = metrics(before)
    after_metrics = metrics(after)
    before_reads = before_metrics["reads"]
    return {
        "before": before_metrics,
        "after": after_metrics,
        "retained_fraction": after_metrics["reads"] / before_reads if before_reads else 0.0,
        "filtering_result": {
            key: int(value)
            for key, value in filtering.items()
            if isinstance(value, (int, float))
        },
        "fastp_version": summary.get(
            "fastp_version", payload.get("fastp_version", "unknown")
        ),
    }


def run_fastp(
    reads: list[Path], output_dir: Path, *, paired: bool, quality: int,
    min_length: int, threads: int,
) -> tuple[list[Path], dict]:
    qc_dir = Path(output_dir) / "preprocessing"
    qc_dir.mkdir(parents=True, exist_ok=True)
    json_path = qc_dir / "fastp.json"
    html_path = qc_dir / "fastp.html"
    cleaned = [qc_dir / "cleaned_R1.fastq.gz"]
    command = [
        "fastp", "--in1", str(reads[0]), "--out1", str(cleaned[0]),
        "--json", str(json_path), "--html", str(html_path),
        "--qualified_quality_phred", str(quality),
        "--length_required", str(min_length), "--thread", str(threads),
    ]
    if paired:
        cleaned.append(qc_dir / "cleaned_R2.fastq.gz")
        command.extend([
            "--in2", str(reads[1]), "--out2", str(cleaned[1]),
            "--detect_adapter_for_pe",
        ])
    returncode, _, stderr = get_formatter().run_subprocess(
        command, "Preprocessing reads with fastp", show_command=True
    )
    if returncode:
        raise MetaQuestError(f"fastp preprocessing failed (exit {returncode}): {stderr}")
    if not json_path.is_file() or not all(path.is_file() for path in cleaned):
        raise MetaQuestError("fastp did not produce its expected cleaned reads and JSON report")
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    result = _summary_metrics(payload)
    result.update({"json_report": str(json_path), "html_report": str(html_path)})
    return cleaned, result
