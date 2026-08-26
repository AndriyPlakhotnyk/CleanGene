from __future__ import annotations

import fcntl
import os
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from .config import truthy

REQUIRED_KRAKEN2_DB_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")
SUPPORTED_KRAKEN2_DATABASE_SIZES = ("standard-8", "standard-16", "standard")


class Kraken2DbError(RuntimeError):
    pass


class Kraken2DbNotReady(Kraken2DbError):
    pass


@dataclass(frozen=True)
class Kraken2DbResolution:
    path: Path
    database_size: str
    source: str
    explicit: bool = False


def normalize_kraken2_database_size(value: str | None) -> str:
    raw = (value or "standard-8").strip().lower()
    aliases = {
        "standard-8": "standard-8",
        "standard_8": "standard-8",
        "8": "standard-8",
        "standard-16": "standard-16",
        "standard_16": "standard-16",
        "16": "standard-16",
        "standard": "standard",
        "full": "standard",
        "full-standard": "standard",
    }
    try:
        return aliases[raw]
    except KeyError:
        supported = "standard-8, standard-16, standard"
        raise Kraken2DbError(f"KRAKEN2_DATABASE_SIZE must be one of {supported} (aliases: standard_8, 8, standard_16, 16, full, full-standard); got {value!r}")


def validate_kraken2_db(path: Path | str) -> None:
    db = Path(path).expanduser()
    if not db.is_dir():
        raise Kraken2DbError(f"Kraken2 database directory not found: {db}")
    for name in REQUIRED_KRAKEN2_DB_FILES:
        required = db / name
        if not required.is_file():
            raise Kraken2DbError(f"Kraken2 database is missing required file: {required}")
        if required.stat().st_size == 0:
            raise Kraken2DbError(f"Kraken2 database required file is empty: {required}")


def kraken2_db_is_valid(path: Path | str) -> bool:
    try:
        validate_kraken2_db(path)
        return True
    except Kraken2DbError:
        return False


def cleangene_project_root() -> Path | None:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pyproject.toml").is_file() and (parent / "src" / "cleangene").is_dir():
            return parent
    return None


def default_kraken2_database_root() -> Path:
    root = cleangene_project_root()
    if root is not None:
        return (root / "databases").resolve()
    return (Path.home() / ".cache" / "cleangene" / "databases").expanduser().resolve()


def kraken2_database_root(cfg: dict[str, str]) -> Path:
    configured = cfg.get("KRAKEN2_DATABASE_ROOT", "").strip()
    return (Path(configured).expanduser() if configured else default_kraken2_database_root()).resolve()


def managed_kraken2_db_path(cfg: dict[str, str]) -> tuple[Path, str, Path]:
    size = normalize_kraken2_database_size(cfg.get("KRAKEN2_DATABASE_SIZE", "standard-8"))
    root = kraken2_database_root(cfg)
    return (root / f"kraken2_{size}").resolve(), size, root


def _script_path() -> Path:
    root = cleangene_project_root()
    if root is None:
        raise Kraken2DbError("Cannot locate scripts/build_kraken2_database.sh; set KRAKEN2_DB to an existing database or run from a source checkout")
    script = root / "scripts" / "build_kraken2_database.sh"
    if not script.is_file():
        raise Kraken2DbError(f"Kraken2 database build script not found: {script}")
    return script


def _publish_staging_database(staging: Path, final: Path) -> None:
    if final.exists():
        if kraken2_db_is_valid(final):
            shutil.rmtree(staging, ignore_errors=True)
            return
        quarantine = final.with_name(f".{final.name}.invalid-{int(time.time())}-{os.getpid()}")
        os.replace(final, quarantine)
    os.replace(staging, final)


def _download_to_shared_database(
    final: Path,
    root: Path,
    size: str,
    cfg: dict[str, str],
    runner: Callable[[Sequence[str]], None],
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    staging = root / f".kraken2_{size}.building"
    if staging.exists() and not staging.is_dir():
        staging.unlink()
    if not staging.exists():
        staging.mkdir(parents=True)
    elif kraken2_db_is_valid(staging):
        shutil.rmtree(staging)
        staging.mkdir(parents=True)
    runner([
        str(_script_path()),
        str(staging),
        cfg.get("KRAKEN2_DB_CPUS", cfg.get("CPUS", "4")),
        cfg.get("KRAKEN2_CLEAN_BUILD_FILES", "true"),
        size,
    ])
    validate_kraken2_db(staging)
    _publish_staging_database(staging, final)
    validate_kraken2_db(final)
    if staging.exists():
        shutil.rmtree(staging, ignore_errors=True)


def resolve_kraken2_db(
    cfg: dict[str, str],
    *,
    allow_download: bool = False,
    runner: Callable[[Sequence[str]], None] | None = None,
    logger: Callable[[str], None] | None = None,
) -> Kraken2DbResolution:
    explicit = cfg.get("KRAKEN2_DB", "").strip()
    size = normalize_kraken2_database_size(cfg.get("KRAKEN2_DATABASE_SIZE", "standard-8"))
    if explicit:
        path = Path(explicit).expanduser().resolve()
        try:
            validate_kraken2_db(path)
        except Kraken2DbError as error:
            raise Kraken2DbError(f"Explicit KRAKEN2_DB is invalid; refusing to auto-download another database: {error}")
        if logger:
            logger(f"Kraken2 database: using explicit KRAKEN2_DB | path={path}")
        return Kraken2DbResolution(path=path, database_size=size, source="explicit", explicit=True)

    path, size, root = managed_kraken2_db_path(cfg)
    try:
        validate_kraken2_db(path)
        if logger:
            logger(f"Kraken2 database: using shared existing database | size={size} | path={path}")
        return Kraken2DbResolution(path=path, database_size=size, source="shared_existing")
    except Kraken2DbError as first_error:
        if not allow_download:
            if not truthy(cfg.get("KRAKEN2_AUTO_DOWNLOAD", "true")):
                raise Kraken2DbNotReady(
                    f"Kraken2 database not found or invalid at {path}: {first_error}. "
                    "KRAKEN2_AUTO_DOWNLOAD=false; create the database there, set KRAKEN2_DATABASE_ROOT, or set a valid KRAKEN2_DB."
                )
            raise Kraken2DbNotReady(
                f"Kraken2 database is not ready at {path}: {first_error}. "
                "The kraken_db_setup stage must complete before preprocess."
            )
        if not truthy(cfg.get("KRAKEN2_AUTO_DOWNLOAD", "true")):
            raise Kraken2DbError(
                f"Kraken2 database not found or invalid at {path}: {first_error}. "
                "KRAKEN2_AUTO_DOWNLOAD=false; create the database there, set KRAKEN2_DATABASE_ROOT, or set a valid KRAKEN2_DB."
            )

    if runner is None:
        raise Kraken2DbError("Internal error: Kraken2 download requested without a runner")
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / f".kraken2_{size}.lock"
    with lock_path.open("w") as lock:
        if logger:
            logger(f"Kraken2 database is being prepared by another CleanGene process; waiting for shared lock | path={lock_path}")
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            validate_kraken2_db(path)
            if logger:
                logger(f"Kraken2 database: using shared existing database | size={size} | path={path}")
            return Kraken2DbResolution(path=path, database_size=size, source="shared_existing")
        except Kraken2DbError:
            if logger:
                logger(f"Kraken2 database not found; preparing shared database | size={size} | path={path}")
            _download_to_shared_database(path, root, size, cfg, runner)
            return Kraken2DbResolution(path=path, database_size=size, source="auto_download")
