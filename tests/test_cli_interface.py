from metaquest.cli import _normalize_global_options, create_parser


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
