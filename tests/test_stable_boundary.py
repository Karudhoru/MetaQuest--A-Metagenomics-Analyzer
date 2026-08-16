import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from metaquest.core import taxonomic_analysis
from metaquest.io import utils
from metaquest.exceptions import AnnotationError, ClassificationError
from metaquest.pipeline.context import PipelineContext
from metaquest.pipeline.runner import PipelineRunner, build_default_pipeline
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.settings import load_config, validate_config


def _stage_names(pipeline):
    return [name for name, _ in pipeline._stages]


def test_experimental_features_are_removed_from_config():
    config = load_config()
    assert not hasattr(config, "ml")
    assert not hasattr(config, "pathogen_detection")


def test_stable_config_validates_without_removed_sections():
    config = load_config()

    assert not hasattr(config, "comparative")
    assert validate_config(config) == (True, [])


def test_legacy_modules_are_removed():
    root = Path(__file__).parents[1] / "src" / "metaquest"
    removed = (
        root / "core" / "pathogen_analysis.py",
        root / "core" / "comparative_analysis.py",
        root / "reporting" / "risk_scoring.py",
        root / "reporting" / "validation_engine.py",
        root / "visualization" / "pathogenic_visualizer.py",
        root / "config" / "pathogen_config.json",
    )

    assert all(not path.exists() for path in removed)


def test_default_pipeline_excludes_experimental_stages():
    config = load_config()

    assert _stage_names(build_default_pipeline(config)) == [
        "Taxonomic Classification",
        "Metagenomic Assembly",
        "Gene Prediction",
        "Functional Annotation",
        "Reporting",
    ]


def test_taxonomy_only_pipeline_skips_assembly_and_annotation():
    config = load_config()

    assert _stage_names(build_default_pipeline(config, skip_annotation=True)) == [
        "Taxonomic Classification",
        "Reporting",
    ]


def test_skip_functional_stops_after_gene_prediction():
    config = load_config()

    assert _stage_names(build_default_pipeline(config, skip_functional=True)) == [
        "Taxonomic Classification",
        "Metagenomic Assembly",
        "Gene Prediction",
        "Reporting",
    ]


def test_failed_functional_stage_preserves_completed_stage_metadata(tmp_path):
    config = load_config(db_dir=tmp_path / "databases")
    context = PipelineContext(config=config, input_files=[], output_dir=tmp_path)
    runner = PipelineRunner(config)

    runner.add_stage("Gene Prediction", lambda ctx: ctx)

    def fail_functional(_ctx):
        raise AnnotationError("emapper failed")

    runner.add_stage("Functional Annotation", fail_functional)

    with pytest.raises(AnnotationError, match="emapper failed"):
        runner.run(context)

    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    assert metadata["status"] == "failed"
    assert metadata["failed_stage"] == "Functional Annotation"
    assert metadata["completed_stages"] == ["Gene Prediction"]


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


def test_low_memory_adds_kraken_memory_mapping(tmp_path, monkeypatch):
    commands = []

    class Formatter:
        def debug(self, _message):
            return None

        def run_subprocess(self, command, **_kwargs):
            commands.append(command)
            (tmp_path / "kraken_report.txt").write_text("report\n", encoding="utf-8")
            return 0, "", ""

    monkeypatch.setattr(taxonomic_analysis, "get_formatter", lambda: Formatter())

    taxonomic_analysis.run_kraken(
        tmp_path / "reads.fastq.gz",
        tmp_path,
        db_path=tmp_path,
        memory_mapping=True,
    )

    assert "--memory-mapping" in commands[0]


def test_gene_prediction_only_check_does_not_require_functional_tools(
    tmp_path, monkeypatch
):
    checked = []
    config = load_config(db_dir=tmp_path)
    config.databases.kraken_db.mkdir(parents=True)
    (config.databases.kraken_db / "hash.k2d").write_text("index\n")

    monkeypatch.setattr(
        utils,
        "_check_command",
        lambda name, _version=None: checked.append(name) or True,
    )

    formatter = SimpleNamespace(
        info=lambda *_args, **_kwargs: None,
        substep=lambda *_args, **_kwargs: None,
        warning=lambda *_args, **_kwargs: None,
        error=lambda *_args, **_kwargs: None,
    )
    utils.run_system_check(formatter, config=config)

    assert "megahit" in checked
    assert "diamond" not in checked
    assert "emapper.py" not in checked


def test_default_functional_check_requires_eggnog_tools_and_databases(
    tmp_path, monkeypatch
):
    checked = []
    config = load_config(db_dir=tmp_path)
    config.databases.kraken_db.mkdir(parents=True)
    (config.databases.kraken_db / "hash.k2d").write_text("index\n")
    config.databases.functional_dir.mkdir(parents=True)
    for name in ("eggnog.db", "eggnog.taxa.db", "eggnog_proteins.dmnd"):
        (config.databases.functional_dir / name).write_text("index\n")

    monkeypatch.setattr(
        utils,
        "_check_command",
        lambda name, _version=None: checked.append(name) or True,
    )
    formatter = SimpleNamespace(
        info=lambda *_args, **_kwargs: None,
        substep=lambda *_args, **_kwargs: None,
        warning=lambda *_args, **_kwargs: None,
        error=lambda *_args, **_kwargs: None,
    )

    utils.run_system_check(formatter, config=config, require_functional=True)

    assert "diamond" in checked
    assert "emapper.py" in checked


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


def test_stable_report_includes_eggnog_provenance(tmp_path):
    functional_dir = tmp_path / "functional_annotation"
    functional_dir.mkdir()
    table = functional_dir / "functional_annotations.tsv"
    categories = functional_dir / "functional_category_summary.tsv"
    table.write_text("query\tannotation_status\ngene_1\tannotated\n")
    categories.write_text("namespace\tterm\tgene_count\nKO\tko:K00001\t1\n")
    (functional_dir / "summary.json").write_text(
        json.dumps(
            {
                "tool_version": "2.1.15",
                "database_release": "5.0.2",
                "tax_scope": "auto",
            }
        )
    )
    context = SimpleNamespace(
        output_dir=tmp_path,
        completed_stages=["Gene Prediction", "Functional Annotation"],
        classification=None,
        assembly=None,
        annotation=SimpleNamespace(
            gene_count=2,
            annotated_count=1,
            functional_annotations=table,
            functional_category_summary=categories,
            functional_reused=True,
        ),
    )

    generate_stable_reports(context)

    summary = json.loads((tmp_path / "analysis_summary.json").read_text())
    assert summary["annotation"]["annotated_genes"] == 1
    assert summary["annotation"]["unannotated_genes"] == 1
    assert summary["annotation"]["eggnog_mapper_version"] == "2.1.15"
    assert summary["annotation"]["eggnog_database_release"] == "5.0.2"
    assert (tmp_path / "03_functional_report.txt").is_file()
