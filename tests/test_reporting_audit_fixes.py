import json
from pathlib import Path
from types import SimpleNamespace

from metaquest.core.preprocessing import _summary_metrics
from metaquest.io.output_formatter import OutputFormatter
from metaquest.pipeline.context import PipelineContext
from metaquest.pipeline.runner import PipelineRunner
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.settings import load_config


def test_fastp_nested_version_and_paired_mean_length():
    result = _summary_metrics(
        {
            "summary": {
                "fastp_version": "0.24.0",
                "before_filtering": {
                    "total_reads": 200,
                    "read1_mean_length": 100,
                    "read2_mean_length": 98,
                },
                "after_filtering": {
                    "total_reads": 150,
                    "read1_mean_length": 95,
                    "read2_mean_length": 93,
                },
            }
        }
    )

    assert result["fastp_version"] == "0.24.0"
    assert result["after"]["mean_length"] == 94


def test_relative_output_and_embedded_html_figures(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "relative-results"
    output.mkdir()
    context = SimpleNamespace(
        output_dir=Path("relative-results"),
        completed_stages=[],
        config=load_config(db_dir=tmp_path / "db"),
        read_mode="single",
        classification=None,
        assembly=None,
        annotation=None,
        preprocessing={
            "before": {"reads": 200, "q20_rate": 0.9, "q30_rate": 0.8},
            "after": {
                "reads": 150,
                "mean_length": 94,
                "q20_rate": 0.98,
                "q30_rate": 0.95,
            },
            "fastp_version": "0.24.0",
            "filtering_result": {"low_quality_reads": 50},
        },
    )

    generate_stable_reports(context)

    report = (output / "report.html").read_text()
    assert "data:image/svg+xml;base64," in report
    assert 'src="plots/' not in report
    assert (output / "analysis_summary.json").is_file()


def test_resume_preserves_log_and_attempt_history(tmp_path):
    log_path = tmp_path / "metaquest.log"
    initial = OutputFormatter(log_file=log_path)
    initial._log("original run")
    initial.log_handle.close()
    resumed = OutputFormatter(log_file=log_path, append_log=True)
    resumed._log("resumed run")
    resumed.log_handle.close()

    log = log_path.read_text()
    assert "original run" in log
    assert "=== MetaQuest resumed run ===" in log
    assert "resumed run" in log

    previous = {
        "status": "completed",
        "started_at": "2026-01-01T00:00:00",
        "timestamp": "2026-01-01T00:10:00",
        "stage_timings_seconds": {"Reporting": 10.0},
        "completed_stages": ["Reporting"],
    }
    (tmp_path / "analysis_metadata.json").write_text(json.dumps(previous))
    config = load_config(db_dir=tmp_path / "databases")
    context = PipelineContext(
        config=config,
        input_files=[],
        output_dir=tmp_path,
        resume=True,
    )
    runner = PipelineRunner(config)
    runner.add_stage("Reporting", lambda ctx: ctx)

    runner.run(context)

    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    assert metadata["started_at"] == "2026-01-01T00:00:00"
    assert metadata["original_stage_timings_seconds"]["Reporting"] == 10.0
    assert metadata["stage_timings_seconds"]["Reporting"] == 10.0
    assert metadata["resume_stage_timings_seconds"]["Reporting"] >= 0
    assert [item["status"] for item in metadata["attempt_history"]] == [
        "completed",
        "completed",
    ]
