"""Deterministic publication-oriented plots and their source TSV files."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def _save(fig, base: Path, formats: tuple[str, ...], dpi: int) -> list[Path]:
    paths = []
    for extension in formats:
        path = base.with_suffix(f".{extension}")
        fig.savefig(path, dpi=dpi if extension == "png" else None, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def _style(palette: str) -> None:
    sns.set_theme(style="whitegrid", palette=palette, context="paper")


def generate_plots(ctx, taxonomy: dict | None) -> list[Path]:
    config = ctx.config.reporting
    output = Path(ctx.output_dir) / "plots"
    data_dir = output / "plot_data"
    data_dir.mkdir(parents=True, exist_ok=True)
    _style(config.color_palette)
    sample = Path(ctx.output_dir).name
    generated: list[Path] = []

    if ctx.preprocessing and ctx.preprocessing.get("before"):
        before, after = ctx.preprocessing["before"], ctx.preprocessing["after"]
        table = pd.DataFrame([
            {"stage": "Before fastp", "metric": "Q20 bases", "percent": before["q20_rate"] * 100},
            {"stage": "After fastp", "metric": "Q20 bases", "percent": after["q20_rate"] * 100},
            {"stage": "Before fastp", "metric": "Q30 bases", "percent": before["q30_rate"] * 100},
            {"stage": "After fastp", "metric": "Q30 bases", "percent": after["q30_rate"] * 100},
        ])
        table.to_csv(data_dir / "qc_metrics.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        sns.barplot(data=table, x="metric", y="percent", hue="stage", ax=ax)
        ax.set(ylim=(0, 100), xlabel="Quality threshold", ylabel="Bases meeting threshold (%)")
        ax.set_title(f"{sample}: read quality before and after preprocessing")
        generated += _save(fig, output / "qc_metrics", config.plot_formats, config.plot_dpi)

    if taxonomy and taxonomy.get("top_taxa"):
        rows = taxonomy["top_taxa"][: config.top_taxa]
        other = max(0.0, 1.0 - sum(row["fraction_all_input_fragments"] for row in rows))
        plot_rows = [{"taxon": row["name"], "percent_all_input_fragments": row["fraction_all_input_fragments"] * 100} for row in rows]
        plot_rows.append({"taxon": "Other/unassigned", "percent_all_input_fragments": other * 100})
        table = pd.DataFrame(plot_rows)
        table.to_csv(data_dir / "taxonomy_top.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(7.2, max(4.2, len(table) * 0.28)))
        sns.barplot(data=table, y="taxon", x="percent_all_input_fragments", color=sns.color_palette()[0], ax=ax)
        ax.set(xlabel="All cleaned input fragments (%)", ylabel="Taxon")
        ax.set_title(f"{sample}: top {min(config.top_taxa, len(rows))} taxa and remainder")
        generated += _save(fig, output / "taxonomy_top", config.plot_formats, config.plot_dpi)

    if ctx.assembly and ctx.assembly.lengths:
        lengths = pd.DataFrame({"contig_length_bp": ctx.assembly.lengths})
        lengths.to_csv(data_dir / "assembly_contig_lengths.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        sns.histplot(data=lengths, x="contig_length_bp", bins="auto", ax=ax)
        ax.set(xlabel="Contig length (bp)", ylabel="Contig count", title=f"{sample}: assembly contig-length distribution")
        generated += _save(fig, output / "assembly_contig_lengths", config.plot_formats, config.plot_dpi)
        cumulative = pd.DataFrame({"contig_length_bp": sorted(ctx.assembly.lengths, reverse=True)})
        cumulative["cumulative_assembly_bp"] = cumulative["contig_length_bp"].cumsum()
        cumulative.insert(0, "contig_rank", range(1, len(cumulative) + 1))
        cumulative.to_csv(data_dir / "assembly_cumulative_length.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        sns.lineplot(data=cumulative, x="contig_rank", y="cumulative_assembly_bp", ax=ax)
        ax.set(xlabel="Contig rank (longest first)", ylabel="Cumulative assembly length (bp)", title=f"{sample}: cumulative assembly length")
        generated += _save(fig, output / "assembly_cumulative_length", config.plot_formats, config.plot_dpi)

    if ctx.annotation and ctx.annotation.functional_category_summary:
        source = pd.read_csv(ctx.annotation.functional_category_summary, sep="\t")
        for namespace in ("COG", "KO", "GO"):
            table = source[source["namespace"] == namespace].nlargest(config.top_functional_terms, "gene_count")
            if table.empty:
                continue
            table.to_csv(data_dir / f"functional_{namespace.lower()}.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(7.2, max(4.2, len(table) * 0.28)))
            sns.barplot(data=table, y="term", x="gene_count", color=sns.color_palette()[0], ax=ax)
            ax.set(xlabel="Annotated genes", ylabel=f"{namespace} term", title=f"{sample}: top {namespace} assignments")
            generated += _save(fig, output / f"functional_{namespace.lower()}", config.plot_formats, config.plot_dpi)
    return generated
