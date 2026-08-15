"""Descriptive reports for the validated MetaQuest pipeline surface."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _native(value: Any) -> Any:
    """Convert pandas/numpy scalar values into JSON-compatible values."""
    if pd.isna(value):
        return None
    return value.item() if hasattr(value, "item") else value


def _taxonomy_summary(bracken_file: Path) -> tuple[dict[str, Any], str]:
    table = pd.read_csv(bracken_file, sep="\t")
    required = {"name", "new_est_reads", "fraction_total_reads"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(
            f"Bracken report is missing required columns: {', '.join(sorted(missing))}"
        )

    table = table.sort_values("fraction_total_reads", ascending=False)
    top = []
    for _, row in table.head(20).iterrows():
        top.append(
            {
                "name": str(row["name"]),
                "estimated_reads": int(row["new_est_reads"]),
                "fraction_total_reads": float(row["fraction_total_reads"]),
            }
        )

    summary = {
        "reported_taxa": int(len(table)),
        "estimated_classified_reads": int(table["new_est_reads"].sum()),
        "top_taxa": top,
    }

    lines = [
        "METAQUEST TAXONOMIC PROFILE",
        "============================",
        "",
        f"Reported taxa: {summary['reported_taxa']}",
        f"Estimated classified reads: {summary['estimated_classified_reads']}",
        "",
        "Top taxa",
        "--------",
        "Taxon\tEstimated reads\tFraction of total reads",
    ]
    lines.extend(
        f"{item['name']}\t{item['estimated_reads']}\t{item['fraction_total_reads']:.6f}"
        for item in top
    )
    lines.extend(
        [
            "",
            "Interpretation note",
            "-------------------",
            "These are computational taxonomic abundance estimates from Kraken2/Bracken.",
            "They are research outputs and do not establish pathogenicity or clinical risk.",
        ]
    )
    return summary, "\n".join(lines) + "\n"


def generate_stable_reports(ctx) -> None:
    """Write descriptive text and JSON outputs without custom risk scoring."""
    output_dir = Path(ctx.output_dir)
    summary: dict[str, Any] = {
        "analysis_type": "short_read_metagenomics",
        "research_use_only": True,
        "completed_stages": list(ctx.completed_stages),
    }

    if ctx.classification:
        taxonomy, report = _taxonomy_summary(ctx.classification.bracken_file)
        summary["taxonomy"] = taxonomy
        (output_dir / "01_taxonomic_report.txt").write_text(report, encoding="utf-8")

    if ctx.assembly:
        summary["assembly"] = {
            "total_contigs": ctx.assembly.total_contigs,
            "total_bases": ctx.assembly.total_bases,
            "n50": ctx.assembly.n50,
            "max_length": ctx.assembly.max_length,
            "mean_length": ctx.assembly.mean_length,
        }

    if ctx.annotation:
        summary["annotation"] = {
            "predicted_genes": ctx.annotation.gene_count,
            "annotated_genes": ctx.annotation.annotated_count,
        }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, default=_native) + "\n",
        encoding="utf-8",
    )
