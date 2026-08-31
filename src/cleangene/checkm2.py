from __future__ import annotations

import fcntl
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from .config import truthy
from .runtime import cleangene_project_root
from .tools import ToolResolutionError, resolve_checkm2_executable
from .util import atomic_json, load_json

EXPECTED_CHECKM2_DB_NAME = "uniref100.KO.1.dmnd"
RUNTIME_VERIFICATION_MARKER = ".cleangene-runtime-verified.json"


class CheckM2DbError(RuntimeError):
    pass


class CheckM2DbNotReady(CheckM2DbError):
    pass


@dataclass(frozen=True)
class CheckM2DbResolution:
    path: Path
    source: str
    explicit: bool = False


def _file_signature(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _checkm2_env_history(executable: Path) -> Path | None:
    for parent in executable.parents:
        if parent.name == "bin":
            candidate = parent.parent / "conda-meta" / "history"
            if candidate.is_file():
                return candidate.resolve()
    return None


def checkm2_runtime_signature(path: Path | str, executable: Path | str, version: str) -> dict[str, object]:
    database = Path(path).expanduser().resolve()
    resolved_executable = Path(executable).expanduser().resolve()
    env_history = _checkm2_env_history(resolved_executable)
    return {
        "schema": 2,
        "CHECKM2_DB": str(database),
        "CHECKM2_DB_FILE": _file_signature(database),
        "CHECKM2_EXECUTABLE": str(resolved_executable),
        "CHECKM2_EXECUTABLE_FILE": _file_signature(resolved_executable),
        "CHECKM2_CONDA_HISTORY_FILE": _file_signature(env_history) if env_history else None,
        "CHECKM2_VERSION": version,
    }


def checkm2_runtime_marker(cfg: dict[str, str]) -> Path:
    return checkm2_database_root(cfg) / RUNTIME_VERIFICATION_MARKER


def checkm2_runtime_is_verified(cfg: dict[str, str], path: Path | str,
                                executable: Path | str, version: str) -> bool:
    marker = checkm2_runtime_marker(cfg)
    if not marker.is_file():
        return False
    try:
        recorded = load_json(marker)
        expected = checkm2_runtime_signature(path, executable, version)
    except (OSError, ValueError):
        return False
    return recorded.get("status") == "complete" and all(recorded.get(key) == value for key, value in expected.items())


def record_checkm2_runtime_verified(cfg: dict[str, str], path: Path | str,
                                    executable: Path | str, version: str) -> Path:
    marker = checkm2_runtime_marker(cfg)
    atomic_json(marker, {"status": "complete", **checkm2_runtime_signature(path, executable, version)})
    return marker


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
        try:
            executable = resolve_checkm2_executable(cfg.get("CHECKM2_EXECUTABLE", ""))
        except ToolResolutionError as error:
            raise CheckM2DbError(f"CheckM2 database download cannot start: {error}") from error
        runner([str(executable), "database", "--download", "--path", str(root), "--no_write_json_db"])
        downloaded = find_managed_checkm2_db(root)
        if not downloaded:
            raise CheckM2DbError(f"CheckM2 download completed but {EXPECTED_CHECKM2_DB_NAME} was not found under {root}")
        validate_checkm2_db(downloaded)
        return CheckM2DbResolution(downloaded, "auto_download")
