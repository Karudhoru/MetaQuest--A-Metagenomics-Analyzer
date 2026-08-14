"""Versioned reference-database installation outside the source tree."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import sqlite3
import subprocess
import tarfile
import tempfile
import urllib.request
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from metaquest.exceptions import MetaQuestError


class DatabaseSetupError(MetaQuestError):
    """A reference database could not be downloaded or validated."""


@dataclass(frozen=True)
class DatabaseSpec:
    key: str
    name: str
    release: str
    status: str
    download_gb: float | None
    installed_gb: float | None
    url: str | None = None
    checksum_url: str | None = None


DATABASES = {
    "taxonomy": DatabaseSpec(
        key="taxonomy",
        name="Kraken2 Standard-8 with Bracken",
        release="2026-06-26",
        status="available",
        download_gb=5.5,
        installed_gb=7.5,
        url=(
            "https://genome-idx.s3.amazonaws.com/kraken/"
            "k2_standard_08_GB_20260626.tar.gz"
        ),
        checksum_url=(
            "https://genome-idx.s3.amazonaws.com/kraken/"
            "standard_08_GB_20260626/standard_08_GB.md5"
        ),
    ),
    "functional": DatabaseSpec(
        key="functional",
        name="eggNOG-mapper core functional database",
        release="5.0.2",
        status="available",
        download_gb=11.3,
        installed_gb=53.0,
        url="http://eggnog5.embl.de/download/emapperdb-5.0.2",
    ),
    "amr": DatabaseSpec(
        "amr", "AMRFinderPlus database", "pending", "planned", None, None,
    ),
    "virulence": DatabaseSpec(
        "virulence", "VFDB Set A", "pending", "planned", None, None,
    ),
}


FUNCTIONAL_ARTIFACTS = (
    {
        "archive": "eggnog.db.gz",
        "output": "eggnog.db",
        "compressed_bytes": 6_776_977_123,
        "format": "gzip",
    },
    {
        "archive": "eggnog.taxa.tar.gz",
        "output": "eggnog.taxa.db",
        "compressed_bytes": 72_797_584,
        "format": "tar",
    },
    {
        "archive": "eggnog_proteins.dmnd.gz",
        "output": "eggnog_proteins.dmnd",
        "compressed_bytes": 5_208_806_170,
        "format": "gzip",
    },
)


def database_rows(db_root: Path) -> list[dict[str, str]]:
    """Return display-ready database status rows."""
    rows = []
    for spec in DATABASES.values():
        installed = is_installed(spec.key, db_root)
        size = (
            f"{spec.download_gb:g} GB download / {spec.installed_gb:g} GB installed"
            if spec.download_gb is not None
            else "not finalized"
        )
        rows.append(
            {
                "Database": spec.key,
                "Description": spec.name,
                "Release": spec.release,
                "Size": size,
                "Status": "installed" if installed else spec.status,
            }
        )
    return rows


def is_installed(key: str, db_root: Path) -> bool:
    """Check the manifest and essential files for an installed database."""
    target = Path(db_root) / key
    manifest = target / "metaquest-db.json"
    if not manifest.exists():
        return False
    if key == "taxonomy":
        return all((target / name).exists() for name in ("hash.k2d", "opts.k2d", "taxo.k2d"))
    if key == "functional":
        return all(
            (target / name).is_file() and (target / name).stat().st_size > 0
            for name in (
                "eggnog.db",
                "eggnog.taxa.db",
                "eggnog_proteins.dmnd",
            )
        )
    return False


def _download(url: str, destination: Path) -> None:
    existing = destination.stat().st_size if destination.exists() else 0
    request = urllib.request.Request(
        url,
        headers={"Range": f"bytes={existing}-"} if existing else {},
    )
    response = urllib.request.urlopen(request)
    resumed = existing and getattr(response, "status", None) == 206
    if resumed:
        content_range = response.headers.get("Content-Range", "")
        if not content_range.startswith(f"bytes {existing}-"):
            response.close()
            raise DatabaseSetupError(
                f"Upstream returned an invalid byte range for {destination.name}"
            )
    mode = "ab" if resumed else "wb"
    with response, destination.open(mode) as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)


def _md5(path: Path) -> str:
    digest = hashlib.md5()  # nosec: upstream publishes MD5 for transfer integrity
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _expected_md5(checksum_file: Path, archive_name: str) -> str:
    """Read the checksum belonging to the archive, not an extracted DB file."""
    for line in checksum_file.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2 and Path(fields[-1]).name == archive_name:
            return fields[0].lower()
    raise DatabaseSetupError(
        f"Upstream checksum file has no entry for {archive_name}"
    )


def _safe_extract(archive: Path, destination: Path) -> None:
    destination_resolved = destination.resolve()
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            member_path = (destination / member.name).resolve()
            if (
                member.issym()
                or member.islnk()
                or destination_resolved not in member_path.parents
            ):
                raise DatabaseSetupError(
                    f"Unsafe path in database archive: {member.name}"
                )
        bundle.extractall(destination)


def _validate_sqlite(path: Path) -> None:
    """Perform a cheap read-only validation without scanning a large database."""
    with path.open("rb") as handle:
        if handle.read(16) != b"SQLite format 3\x00":
            raise DatabaseSetupError(f"Invalid SQLite database header: {path.name}")
    try:
        connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            connection.execute("SELECT name FROM sqlite_master LIMIT 1").fetchone()
        finally:
            connection.close()
    except sqlite3.DatabaseError as exc:
        raise DatabaseSetupError(f"Invalid SQLite database: {path.name}") from exc


def _validate_functional_tools(target: Path, release: str) -> None:
    diamond = shutil.which("diamond")
    emapper = shutil.which("emapper.py")
    if not diamond or not emapper:
        raise DatabaseSetupError(
            "DIAMOND and eggNOG-mapper are required to validate functional data"
        )

    commands = (
        ([diamond, "dbinfo", "--db", str(target / "eggnog_proteins.dmnd")], "DIAMOND"),
        ([emapper, "--data_dir", str(target), "--version"], "eggNOG-mapper"),
    )
    for command, name in commands:
        result = subprocess.run(command, capture_output=True, text=True)
        output = f"{result.stdout}\n{result.stderr}"
        if result.returncode or (name == "eggNOG-mapper" and release not in output):
            raise DatabaseSetupError(f"{name} could not validate the functional database")


def _extract_functional_artifact(
    archive: Path,
    destination: Path,
    artifact: dict[str, object],
) -> None:
    if artifact["format"] == "tar":
        _safe_extract(archive, destination)
        return
    output = destination / str(artifact["output"])
    try:
        with gzip.open(archive, "rb") as source, output.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    except (EOFError, OSError) as exc:
        raise DatabaseSetupError(f"Invalid gzip archive: {archive.name}") from exc


def _install_functional(
    spec: DatabaseSpec,
    db_root: Path,
    target: Path,
    *,
    force: bool,
    notify: Callable[[str], None] | None,
) -> Path:
    downloads = db_root / ".downloads" / "functional"
    downloads.mkdir(parents=True, exist_ok=True)
    provenance = []

    for artifact in FUNCTIONAL_ARTIFACTS:
        archive_name = str(artifact["archive"])
        partial = downloads / f"{archive_name}.partial"
        url = f"{spec.url}/{archive_name}"
        if notify:
            if partial.exists() and partial.stat().st_size:
                notify(
                    f"Downloading {archive_name} from existing partial file "
                    f"({partial.stat().st_size / 1024**3:.1f} GB)"
                )
            else:
                notify(f"Downloading {archive_name}")
        expected_bytes = int(artifact["compressed_bytes"])
        if not partial.exists() or partial.stat().st_size != expected_bytes:
            _download(url, partial)
        if partial.stat().st_size != expected_bytes:
            raise DatabaseSetupError(
                f"Unexpected size for {archive_name}: expected {expected_bytes} bytes, "
                f"got {partial.stat().st_size}"
            )
        if notify:
            notify(f"Calculating SHA-256 for {archive_name}")
        provenance.append(
            {
                "url": url,
                "compressed_bytes": expected_bytes,
                "sha256": _sha256(partial),
            }
        )

    with tempfile.TemporaryDirectory(prefix=".metaquest-functional-", dir=db_root) as tmp:
        extracted = Path(tmp) / "functional"
        extracted.mkdir()
        for artifact in FUNCTIONAL_ARTIFACTS:
            archive_name = str(artifact["archive"])
            if notify:
                notify(f"Extracting {archive_name}")
            _extract_functional_artifact(
                downloads / f"{archive_name}.partial",
                extracted,
                artifact,
            )

        required = tuple(str(item["output"]) for item in FUNCTIONAL_ARTIFACTS)
        for name in required:
            path = extracted / name
            if not path.is_file() or path.stat().st_size == 0:
                raise DatabaseSetupError(f"Functional database is missing required file: {name}")
        _validate_sqlite(extracted / "eggnog.db")
        _validate_sqlite(extracted / "eggnog.taxa.db")
        _validate_functional_tools(extracted, spec.release)

        manifest = {
            **asdict(spec),
            "installed_at": datetime.now(timezone.utc).isoformat(),
            "artifacts": provenance,
            "source": "eggNOG v5 official download server",
        }
        (extracted / "metaquest-db.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

        previous = None
        if target.exists():
            previous = Path(tmp) / "previous"
            target.replace(previous)
        try:
            extracted.replace(target)
        except Exception:
            if previous is not None and previous.exists():
                previous.replace(target)
            raise

    for artifact in FUNCTIONAL_ARTIFACTS:
        (downloads / f"{artifact['archive']}.partial").unlink()
    try:
        downloads.rmdir()
        downloads.parent.rmdir()
    except OSError:
        pass
    if notify:
        notify("Functional database manifest written and installation validated")
    return target


def install_database(
    key: str,
    db_root: Path,
    *,
    force: bool = False,
    notify: Callable[[str], None] | None = None,
) -> Path:
    """Download, verify, extract, validate, and atomically install a database."""
    if key not in DATABASES:
        raise DatabaseSetupError(f"Unknown database: {key}")
    spec = DATABASES[key]
    if spec.status != "available" or not spec.url:
        raise DatabaseSetupError(
            f"Database '{key}' is planned but its installation is not finalized"
        )

    db_root = Path(db_root).expanduser().resolve()
    target = db_root / key
    if is_installed(key, db_root) and not force:
        if notify:
            notify(f"{spec.name} is already installed and valid")
        return target
    if target.exists() and not force:
        raise DatabaseSetupError(
            f"Database directory exists but is not valid: {target}",
            suggestions=["Move it aside or rerun with --force to replace it."],
        )

    required_bytes = int((spec.download_gb + spec.installed_gb + 1) * 1024**3)
    db_root.mkdir(parents=True, exist_ok=True)
    if notify:
        notify(
            f"Disk preflight: {spec.download_gb:g} GB download, "
            f"{spec.installed_gb:g} GB installed"
        )
    if shutil.disk_usage(db_root).free < required_bytes:
        raise DatabaseSetupError(
            f"Insufficient free space for {spec.name}; "
            f"at least {required_bytes / 1024**3:.1f} GB is required"
        )

    if key == "functional":
        return _install_functional(
            spec,
            db_root,
            target,
            force=force,
            notify=notify,
        )
    if not spec.checksum_url:
        raise DatabaseSetupError(f"No checksum source is configured for {spec.name}")

    downloads = db_root / ".downloads"
    downloads.mkdir(exist_ok=True)
    archive_name = Path(spec.url).name
    archive = downloads / f"{archive_name}.partial"
    checksum_file = downloads / f"{archive_name}.md5"

    if notify and archive.exists():
        notify(
            f"Resuming database archive at "
            f"{archive.stat().st_size / 1024**3:.1f} GB"
        )

    if notify:
        notify("Downloading upstream checksum")
    checksum_file.unlink(missing_ok=True)
    _download(spec.checksum_url, checksum_file)
    if notify:
        notify("Downloading database archive")
    _download(spec.url, archive)
    if notify:
        notify("Verifying archive checksum")
    expected = _expected_md5(checksum_file, archive_name)
    observed = _md5(archive)
    if observed != expected:
        archive.unlink(missing_ok=True)
        raise DatabaseSetupError(
            f"Checksum mismatch for {spec.name}: expected {expected}, got {observed}"
        )

    with tempfile.TemporaryDirectory(prefix=".metaquest-db-", dir=db_root) as tmp:
        work = Path(tmp)
        extracted = work / "extracted"
        extracted.mkdir()

        if notify:
            notify("Extracting database archive")
        _safe_extract(archive, extracted)
        required = ("hash.k2d", "opts.k2d", "taxo.k2d")
        if not all((extracted / name).exists() for name in required):
            raise DatabaseSetupError(
                f"{spec.name} archive is missing required Kraken2 index files"
            )

        manifest = {
            **asdict(spec),
            "installed_at": datetime.now(timezone.utc).isoformat(),
            "archive_md5": observed,
            "source": "https://benlangmead.github.io/aws-indexes/k2",
        }
        (extracted / "metaquest-db.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

        if target.exists():
            shutil.rmtree(target)
        extracted.replace(target)
    archive.unlink(missing_ok=True)
    checksum_file.unlink(missing_ok=True)
    try:
        downloads.rmdir()
    except OSError:
        pass
    if notify:
        notify("Database manifest written and installation validated")
    return target
