from metaquest.cli import _normalize_global_options, _prepare_output_directory, create_parser
import pytest


def test_modern_run_command_and_legacy_analyze_alias_match():
    parser = create_parser()
    modern = parser.parse_args(['run', '--single', 'reads.fastq.gz'])
    legacy = parser.parse_args(['analyze', '--single', 'reads.fastq.gz'])

    assert modern.command == legacy.command == 'run'
    assert modern.single == legacy.single == 'reads.fastq.gz'


def test_database_command_and_legacy_alias_match():
    parser = create_parser()

    assert parser.parse_args(['databases', '--list']).command == 'databases'
    assert parser.parse_args(['setup-db', '--list']).command == 'databases'


def test_compare_stub_does_not_require_legacy_inputs():
    assert create_parser().parse_args(['compare']).command == 'compare'


def test_global_output_options_can_follow_command():
    normalized = _normalize_global_options(
        ['run', '--quiet', '--single', 'reads.fastq.gz', '--no-color']
    )
    args = create_parser().parse_args(normalized)

    assert args.command == 'run'
    assert args.quiet is True
    assert args.no_color is True


def test_low_memory_flag_can_precede_or_follow_run_command():
    parser = create_parser()

    before = parser.parse_args(_normalize_global_options(
        ['--low-memory', 'run', '--single', 'reads.fastq.gz']
    ))
    after = parser.parse_args(_normalize_global_options(
        ['run', '--low-memory', '--single', 'reads.fastq.gz']
    ))

    assert before.low_memory is True
    assert after.low_memory is True


def test_functional_workflow_controls_are_explicit():
    parser = create_parser()

    default = parser.parse_args(["run", "--single", "reads.fastq.gz"])
    taxonomy = parser.parse_args(
        ["run", "--taxonomy-only", "--single", "reads.fastq.gz"]
    )
    genes_only = parser.parse_args(
        ["run", "--skip-functional", "--single", "reads.fastq.gz"]
    )

    assert default.taxonomy_only is False and default.skip_functional is False
    assert taxonomy.taxonomy_only is True
    assert genes_only.skip_functional is True


def test_force_and_resume_are_explicit_and_mutually_exclusive():
    parser = create_parser()

    assert parser.parse_args(["run", "--force", "--single", "reads.fastq"]).force
    assert parser.parse_args(["run", "--resume", "--single", "reads.fastq"]).resume
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["run", "--force", "--resume", "--single", "reads.fastq"]
        )


def test_force_preserves_existing_output_as_backup(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    (output / "analysis_metadata.json").write_text("{}")
    parser = create_parser()
    args = parser.parse_args(
        ["run", "--force", "--output", str(output), "--single", "reads.fastq"]
    )

    _prepare_output_directory(args, parser)

    assert output.is_dir() and not any(output.iterdir())
    backups = list(tmp_path.glob("results.metaquest-backup-*"))
    assert len(backups) == 1
    assert (backups[0] / "analysis_metadata.json").is_file()
