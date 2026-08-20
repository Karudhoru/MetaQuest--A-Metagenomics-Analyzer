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


def _taxonomy_summary(classification) -> tuple[dict[str, Any], str]:
    bracken_file = classification.bracken_file
    table = pd.read_csv(bracken_file, sep="\t")
    required = {"name", "new_est_reads"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(
            f"Bracken report is missing required columns: {', '.join(sorted(missing))}"
        )

    species_assigned = classification.species_assigned_reads
    total_reads = classification.total_reads
    table = table.sort_values("new_est_reads", ascending=False)
    top = []
    for _, row in table.head(20).iterrows():
        top.append(
            {
                "name": str(row["name"]),
                "estimated_reads": int(row["new_est_reads"]),
                "fraction_species_assigned_reads": (
                    int(row["new_est_reads"]) / species_assigned if species_assigned else 0.0
                ),
                "fraction_all_input_reads": (
                    int(row["new_est_reads"]) / total_reads if total_reads else 0.0
                ),
            }
        )

    summary = {
        "reported_taxa": int(len(table)),
        "total_input_reads": total_reads,
        "kraken_classified_reads": classification.kraken_classified_reads,
        "unclassified_reads": classification.unclassified_reads,
        "classification_rate": (
            classification.kraken_classified_reads / total_reads if total_reads else 0.0
        ),
        "species_assigned_reads": species_assigned,
        "observed_read_length": classification.observed_read_length,
        "bracken_read_length": classification.bracken_read_length,
        "top_taxa": top,
    }

    lines = [
        "METAQUEST TAXONOMIC PROFILE",
        "============================",
        "",
        f"Reported taxa: {summary['reported_taxa']}",
        f"Total input reads: {summary['total_input_reads']}",
        f"Kraken-classified reads: {summary['kraken_classified_reads']}",
        f"Unclassified reads: {summary['unclassified_reads']}",
        f"Classification rate: {summary['classification_rate']:.2%}",
        f"Species-assigned reads: {summary['species_assigned_reads']}",
        f"Observed read length: {summary['observed_read_length']} bp",
        f"Bracken model length: {summary['bracken_read_length']} bp",
        "",
        "Top taxa",
        "--------",
        "Taxon\tEstimated reads\tFraction species-assigned\tFraction all input",
    ]
    lines.extend(
        f"{item['name']}\t{item['estimated_reads']}\t"
        f"{item['fraction_species_assigned_reads']:.6f}\t"
        f"{item['fraction_all_input_reads']:.6f}"
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
    output_dir = Path(ctx.output_dir).resolve()
    summary: dict[str, Any] = {
        "analysis_type": "short_read_metagenomics",
        "research_use_only": True,
        "completed_stages": [*ctx.completed_stages, "Reporting"],
    }

    if ctx.classification:
        taxonomy, report = _taxonomy_summary(ctx.classification)
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
        assembly_lines = [
            "METAQUEST METAGENOMIC ASSEMBLY",
            "===============================",
            "",
            f"Contigs: {ctx.assembly.total_contigs}",
            f"Total bases: {ctx.assembly.total_bases}",
            f"N50: {ctx.assembly.n50} bp",
            f"Maximum contig length: {ctx.assembly.max_length} bp",
            f"Mean contig length: {ctx.assembly.mean_length:.2f} bp",
            "",
            "Interpretation note",
            "-------------------",
            "Assembly continuity does not by itself establish assembly correctness.",
        ]
        (output_dir / "02_assembly_report.txt").write_text(
            "\n".join(assembly_lines) + "\n", encoding="utf-8"
        )

    if ctx.annotation:
        annotation_summary = {
            "predicted_genes": ctx.annotation.gene_count,
            "annotated_genes": ctx.annotation.annotated_count,
            "unannotated_genes": ctx.annotation.gene_count - ctx.annotation.annotated_count,
            "functional_status": (
                "completed" if ctx.annotation.functional_annotations else "skipped"
            ),
        }
        if ctx.annotation.functional_annotations:
            functional_dir = ctx.annotation.functional_annotations.parent
            eggnog_summary_path = functional_dir / "summary.json"
            eggnog_summary = json.loads(eggnog_summary_path.read_text(encoding="utf-8"))
            annotation_summary.update(
                {
                    "functional_annotations": str(
                        ctx.annotation.functional_annotations.relative_to(output_dir)
                    ),
                    "functional_category_summary": str(
                        ctx.annotation.functional_category_summary.relative_to(output_dir)
                    ),
                    "functional_reused": ctx.annotation.functional_reused,
                    "eggnog_mapper_version": eggnog_summary["tool_version"],
                    "eggnog_database_release": eggnog_summary["database_release"],
                    "tax_scope": eggnog_summary["tax_scope"],
                }
            )
            lines = [
                "METAQUEST FUNCTIONAL ANNOTATION",
                "===============================",
                "",
                f"Predicted genes: {ctx.annotation.gene_count}",
                f"eggNOG-annotated genes: {ctx.annotation.annotated_count}",
                f"Unannotated genes: {ctx.annotation.gene_count - ctx.annotation.annotated_count}",
                f"eggNOG-mapper: {eggnog_summary['tool_version']}",
                f"eggNOG database: {eggnog_summary['database_release']}",
                f"Taxonomic scope: {eggnog_summary['tax_scope']}",
                "",
                "Interpretation note",
                "-------------------",
                "Annotations are computational orthology-based assignments for research use.",
                "They do not establish gene expression, phenotype, or clinical significance.",
            ]
            (output_dir / "03_functional_report.txt").write_text(
                "\n".join(lines) + "\n", encoding="utf-8"
            )
        summary["annotation"] = annotation_summary
    summary_path = output_dir / "analysis_summary.json"
    temporary_summary = summary_path.with_suffix(".json.tmp")
    temporary_summary.write_text(
        json.dumps(summary, indent=2, default=_native) + "\n",
        encoding="utf-8",
    )
    temporary_summary.replace(summary_path)
