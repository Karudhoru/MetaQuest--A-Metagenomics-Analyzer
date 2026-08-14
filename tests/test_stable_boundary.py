from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from metaquest.core import taxonomic_analysis
from metaquest.exceptions import ClassificationError
from metaquest.pipeline.runner import build_default_pipeline
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.settings import load_config


def _stage_names(pipeline):
    return [name for name, _ in pipeline._stages]


def test_experimental_features_are_removed_from_config():
    config = load_config()
    assert not hasattr(config, "ml")
    assert not hasattr(config, "pathogen_detection")


def test_default_pipeline_excludes_experimental_stages():
    config = load_config()

    assert _stage_names(build_default_pipeline(config)) == [
        "Taxonomic Classification",
        "Metagenomic Assembly",
        "Functional Annotation",
        "Reporting",
    ]


def test_taxonomy_only_pipeline_skips_assembly_and_annotation():
    config = load_config()

    assert _stage_names(build_default_pipeline(config, skip_annotation=True)) == [
        "Taxonomic Classification",
        "Reporting",
    ]


def test_bracken_failure_is_not_treated_as_a_bracken_report(tmp_path, monkeypatch):
    kraken_report = tmp_path / "kraken_report.txt"
    kraken_report.write_text("report\n", encoding="utf-8")

    class Formatter:
        def debug(self, _message):
            return None

        def run_subprocess(self, *_args, **_kwargs):
            return 2, "", "bracken failed"

    monkeypatch.setattr(taxonomic_analysis, "get_formatter", lambda: Formatter())

    with pytest.raises(ClassificationError, match="Bracken failed"):
        taxonomic_analysis.run_bracken(
            kraken_report,
            tmp_path,
            db_path=tmp_path,
        )


def test_stable_report_contains_no_risk_output(tmp_path):
    bracken_file = tmp_path / "bracken_report.tsv"
    pd.DataFrame(
        {
            "name": ["Taxon A", "Taxon B"],
            "new_est_reads": [80, 20],
            "fraction_total_reads": [0.8, 0.2],
        }
    ).to_csv(bracken_file, sep="\t", index=False)

    context = SimpleNamespace(
        output_dir=tmp_path,
        completed_stages=["Taxonomic Classification"],
        classification=SimpleNamespace(bracken_file=bracken_file),
        assembly=None,
        annotation=None,
    )

    generate_stable_reports(context)

    summary = (tmp_path / "analysis_summary.json").read_text(encoding="utf-8")
    report = (tmp_path / "01_taxonomic_report.txt").read_text(encoding="utf-8")

    assert "risk_score" not in summary
    assert "clinical recommendation" not in report.lower()
    assert "do not establish pathogenicity or clinical risk" in report
