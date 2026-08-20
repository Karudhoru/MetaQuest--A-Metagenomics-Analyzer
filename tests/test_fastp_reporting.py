from metaquest.core.preprocessing import _summary_metrics
from metaquest.pipeline.stages.assembly import _compute_assembly_stats
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.settings import load_config
from types import SimpleNamespace


def test_fastp_summary_uses_rates_and_retention():
    result = _summary_metrics({
        "fastp_version": "0.24.0",
        "summary": {
            "before_filtering": {"total_reads": 200, "total_bases": 20000, "q20_rate": 0.9, "q30_rate": 0.8, "gc_content": 0.5},
            "after_filtering": {"total_reads": 150, "total_bases": 14000, "q20_rate": 0.98, "q30_rate": 0.95, "gc_content": 0.49},
        },
        "filtering_result": {"low_quality_reads": 40, "too_short_reads": 10},
    })
    assert result["retained_fraction"] == 0.75
    assert result["after"]["q30_rate"] == 0.95
    assert result["filtering_result"]["low_quality_reads"] == 40


def test_expanded_assembly_statistics(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">a\n" + "G" * 1000 + "\n>b\n" + "A" * 500 + "\n>c\n" + "C" * 100 + "\n")
    stats = _compute_assembly_stats(fasta)
    assert stats["n50"] == 1000
    assert stats["l50"] == 1
    assert stats["n90"] == 500
    assert stats["l90"] == 2
    assert stats["contigs_ge_1000"] == 1
    assert stats["gc_percent"] == 68.75


def test_offline_report_and_plot_data_are_generated(tmp_path):
    ctx = SimpleNamespace(
        output_dir=tmp_path,
        completed_stages=[],
        config=load_config(db_dir=tmp_path / "db"),
        read_mode="single",
        classification=None,
        assembly=None,
        annotation=None,
        preprocessing={
            "before": {"q20_rate": 0.90, "q30_rate": 0.80},
            "after": {"q20_rate": 0.99, "q30_rate": 0.95},
            "retained_fraction": 0.8,
        },
    )
    generate_stable_reports(ctx)
    assert (tmp_path / "report.html").is_file()
    assert (tmp_path / "plots" / "qc_metrics.svg").is_file()
    assert (tmp_path / "plots" / "qc_metrics.png").is_file()
    assert (tmp_path / "plots" / "plot_data" / "qc_metrics.tsv").is_file()
