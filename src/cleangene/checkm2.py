from __future__ import annotations

import fcntl
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from .config import truthy
from .runtime import cleangene_project_root

EXPECTED_CHECKM2_DB_NAME = "uniref100.KO.1.dmnd"


class CheckM2DbError(RuntimeError):
    pass


class CheckM2DbNotReady(CheckM2DbError):
    pass


@dataclass(frozen=True)
class CheckM2DbResolution:
    path: Path
    source: str
    explicit: bool = False


def validate_checkm2_db(path: Path | str) -> None:
    db = Path(path).expanduser()
    if not db.is_file():
        raise CheckM2DbError(f"CheckM2 database file not found: {db}")
    if db.stat().st_size == 0:
        raise CheckM2DbError(f"CheckM2 database file is empty: {db}")
    if db.name != EXPECTED_CHECKM2_DB_NAME:
        raise CheckM2DbError(f"CheckM2 database must be named {EXPECTED_CHECKM2_DB_NAME}: {db}")


def checkm2_db_is_valid(path: Path | str) -> bool:
    try:
        validate_checkm2_db(path)
        return True
    except CheckM2DbError:
        return False


def default_checkm2_database_root() -> Path:
    root = cleangene_project_root()
    if root is not None:
        return (root / "databases" / "checkm2").resolve()
    return (Path.home() / ".cache" / "cleangene" / "databases" / "checkm2").expanduser().resolve()


def checkm2_database_root(cfg: dict[str, str]) -> Path:
    configured = cfg.get("CHECKM2_DATABASE_ROOT", "").strip()
    generic = cfg.get("CLEANGENE_DATABASE_ROOT", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    if generic:
        return (Path(generic).expanduser() / "checkm2").resolve()
    return default_checkm2_database_root()


def find_managed_checkm2_db(root: Path) -> Path | None:
    direct = root / "CheckM2_database" / EXPECTED_CHECKM2_DB_NAME
    if checkm2_db_is_valid(direct):
        return direct.resolve()
    flat = root / EXPECTED_CHECKM2_DB_NAME
    if checkm2_db_is_valid(flat):
        return flat.resolve()
    candidates = [p.resolve() for p in root.rglob(EXPECTED_CHECKM2_DB_NAME) if checkm2_db_is_valid(p)] if root.is_dir() else []
    unique = sorted(set(candidates))
    if len(unique) == 1:
        return unique[0]
    if len(unique) > 1:
        raise CheckM2DbError("Multiple CheckM2 database candidates found under " + str(root) + ": " + ", ".join(map(str, unique[:5])))
    return None


def resolve_checkm2_db(
    cfg: dict[str, str],
    *,
    allow_download: bool = False,
    runner: Callable[[Sequence[str]], None] | None = None,
    logger: Callable[[str], None] | None = None,
) -> CheckM2DbResolution:
    explicit = cfg.get("CHECKM2_DB", "").strip()
    if explicit:
        path = Path(explicit).expanduser().resolve()
        try:
            validate_checkm2_db(path)
        except CheckM2DbError as error:
            raise CheckM2DbError(f"Explicit CHECKM2_DB is invalid; refusing to auto-download another database: {error}")
        if logger:
            logger(f"CheckM2 database: using explicit CHECKM2_DB | path={path}")
        return CheckM2DbResolution(path, "explicit", True)

    root = checkm2_database_root(cfg)
    try:
        existing = find_managed_checkm2_db(root)
    except CheckM2DbError:
        raise
    if existing:
        if logger:
            logger(f"CheckM2 database: using shared existing database | path={existing}")
        return CheckM2DbResolution(existing, "shared_existing")

    if not allow_download:
        message = (
            f"CheckM2 database is not ready under {root}. "
            "The checkm2_db_setup stage must complete before preprocess."
        )
        if not truthy(cfg.get("CHECKM2_AUTO_DOWNLOAD", "true")):
            message = (
                f"CheckM2 database not found under {root}. CHECKM2_AUTO_DOWNLOAD=false; "
                "create the database there, set CHECKM2_DATABASE_ROOT, or set a valid CHECKM2_DB."
            )
        raise CheckM2DbNotReady(message)
    if not truthy(cfg.get("CHECKM2_AUTO_DOWNLOAD", "true")):
        raise CheckM2DbError(
            f"CheckM2 database not found under {root}. CHECKM2_AUTO_DOWNLOAD=false; "
            "create the database there, set CHECKM2_DATABASE_ROOT, or set a valid CHECKM2_DB."
        )
    if runner is None:
        raise CheckM2DbError("Internal error: CheckM2 download requested without a runner")

    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / ".download.lock"
    with lock_path.open("w") as lock:
        if logger:
            logger(f"CheckM2 database is being prepared by another CleanGene process; waiting for shared lock | path={lock_path}")
        fcntl.flock(lock, fcntl.LOCK_EX)
        existing = find_managed_checkm2_db(root)
        if existing:
            if logger:
                logger(f"CheckM2 database: using shared existing database | path={existing}")
            return CheckM2DbResolution(existing, "shared_existing")
        if logger:
            logger(f"CheckM2 database not found; preparing shared database | path={root}")
        executable = cfg.get("CHECKM2_EXECUTABLE", "").strip() or "checkm2"
        runner([executable, "database", "--download", "--path", str(root), "--no_write_json_db"])
        downloaded = find_managed_checkm2_db(root)
        if not downloaded:
            raise CheckM2DbError(f"CheckM2 download completed but {EXPECTED_CHECKM2_DB_NAME} was not found under {root}")
        validate_checkm2_db(downloaded)
        return CheckM2DbResolution(downloaded, "auto_download")
