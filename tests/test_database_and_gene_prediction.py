import io
import gzip
import sqlite3
import sys
import tarfile
from types import SimpleNamespace

import pytest

from metaquest.core.functional_analysis import run_gene_prediction
from metaquest.database_manager import (
    DatabaseSetupError,
    _extract_functional_artifact,
    _expected_md5,
    _safe_extract,
    database_rows,
)


def test_database_catalog_distinguishes_available_and_planned(tmp_path):
    rows = {row["Database"]: row for row in database_rows(tmp_path)}

    assert rows["taxonomy"]["Status"] == "available"
    assert rows["taxonomy"]["Release"] == "2026-06-26"
    assert rows["functional"]["Status"] == "available"
    assert rows["functional"]["Release"] == "5.0.2"
    assert rows["amr"]["Status"] == "planned"
    assert rows["virulence"]["Status"] == "planned"


def test_safe_extract_rejects_path_traversal(tmp_path):
    archive = tmp_path / "unsafe.tar.gz"
    payload = b"not allowed"
    with tarfile.open(archive, "w:gz") as bundle:
        member = tarfile.TarInfo("../outside.txt")
        member.size = len(payload)
        bundle.addfile(member, io.BytesIO(payload))

    destination = tmp_path / "extract"
    destination.mkdir()
    with pytest.raises(DatabaseSetupError, match="Unsafe path"):
        _safe_extract(archive, destination)

    assert not (tmp_path / "outside.txt").exists()


def test_archive_checksum_is_selected_by_filename(tmp_path):
    checksums = tmp_path / "checksums.md5"
    checksums.write_text(
        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa  hash.k2d\n"
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb  database.tar.gz\n",
        encoding="utf-8",
    )

    assert _expected_md5(checksums, "database.tar.gz") == "b" * 32


def test_functional_gzip_extraction_streams_to_expected_name(tmp_path):
    archive = tmp_path / "eggnog.db.gz.partial"
    with gzip.open(archive, "wb") as handle:
        handle.write(b"SQLite format 3\x00test payload")
    destination = tmp_path / "output"
    destination.mkdir()

    _extract_functional_artifact(
        archive,
        destination,
        {"format": "gzip", "output": "eggnog.db"},
    )

    assert (destination / "eggnog.db").read_bytes() == b"SQLite format 3\x00test payload"


def test_functional_taxonomy_archive_uses_safe_extraction(tmp_path):
    source = tmp_path / "source.db"
    with sqlite3.connect(source) as connection:
        connection.execute("CREATE TABLE taxa (id INTEGER)")
    archive = tmp_path / "eggnog.taxa.tar.gz.partial"
    with tarfile.open(archive, "w:gz") as bundle:
        bundle.add(source, arcname="eggnog.taxa.db")
    destination = tmp_path / "output"
    destination.mkdir()

    _extract_functional_artifact(
        archive,
        destination,
        {"format": "tar", "output": "eggnog.taxa.db"},
    )

    assert (destination / "eggnog.taxa.db").is_file()


def test_pyrodigal_outputs_are_streamed_per_contig(tmp_path, monkeypatch):
    class FakePredictions:
        def __len__(self):
            return 1

        def write_translations(self, handle, sequence_id, **_kwargs):
            handle.write(f">{sequence_id}_1\nMPEPTIDE\n")

        def write_genes(self, handle, sequence_id, **_kwargs):
            handle.write(f">{sequence_id}_1\nATGCCCTAA\n")

        def write_gff(self, handle, sequence_id, **_kwargs):
            handle.write(
                f"{sequence_id}\tPyrodigal\tCDS\t1\t9\t.\t+\t0\tID={sequence_id}_1\n"
            )

    class FakeGeneFinder:
        def __init__(self, *, meta):
            assert meta is True

        def find_genes(self, sequence):
            assert isinstance(sequence, bytes)
            return FakePredictions()

    monkeypatch.setitem(
        sys.modules,
        "pyrodigal",
        SimpleNamespace(GeneFinder=FakeGeneFinder, __version__="test"),
    )

    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(
        ">long_contig\n" + "ATG" * 100 + "\n>short_contig\nATG\n",
        encoding="utf-8",
    )

    output_dir, count = run_gene_prediction(
        contigs,
        tmp_path / "results",
        min_contig_length=200,
    )

    assert count == 1
    assert ">long_contig_1" in (output_dir / "genes.faa").read_text()
    assert "short_contig" not in (output_dir / "genes.faa").read_text()
    assert (output_dir / "genes.fna").exists()
    assert (output_dir / "genes.gff3").read_text().startswith("##gff-version 3")
    assert (output_dir / "summary.json").exists()
