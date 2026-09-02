from __future__ import annotations

import fcntl
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from .config import truthy
from .runtime import cleangene_project_root
from .tools import ToolResolutionError, resolve_checkm2_executable
from .util import atomic_json, load_json, read_tsv, safe_name

EXPECTED_CHECKM2_DB_NAME = "uniref100.KO.1.dmnd"
RUNTIME_VERIFICATION_MARKER = ".cleangene-runtime-verified.json"
CHECKM2_COMMAND_SCHEMA_VERSION = 4


class CheckM2DbError(RuntimeError):
    pass


class CheckM2DbNotReady(CheckM2DbError):
    pass


@dataclass(frozen=True)
class CheckM2DbResolution:
    path: Path
    source: str
    explicit: bool = False


@dataclass(frozen=True)
class CheckM2PredictCapabilities:
    cleanup_option: str
    help_sha256: str


def _sha256_text(text: str) -> str:
    import hashlib
    return hashlib.sha256(text.encode()).hexdigest()


def checkm2_predict_help(executable: Path | str, timeout_seconds: int = 60) -> str:
    command = [str(executable), "predict", "--help"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout_seconds)
    except subprocess.TimeoutExpired as error:
        raise CheckM2DbError(f"CheckM2 predict help timed out after {timeout_seconds} seconds") from error
    except (OSError, subprocess.SubprocessError) as error:
        raise CheckM2DbError(f"CheckM2 predict help could not run: {error}") from error
    text = (result.stdout or "") + "\n" + (result.stderr or "")
    if result.returncode:
        raise CheckM2DbError(f"CheckM2 predict help failed with exit status {result.returncode}: {text.strip()}")
    return text


def checkm2_predict_capabilities(executable: Path | str, help_text: str | None = None) -> CheckM2PredictCapabilities:
    help_text = help_text if help_text is not None else checkm2_predict_help(executable)
    required = ("--input", "--output-directory", "--database_path", "--threads", "--force")
    missing = [option for option in required if option not in help_text]
    if missing:
        raise CheckM2DbError("CheckM2 predict CLI is missing required option(s): " + ", ".join(missing))
    if "--remove_intermediates" in help_text:
        cleanup = "--remove_intermediates"
    else:
        future_cleanup = "--remove" + "-intermediates"
        if future_cleanup in help_text:
            cleanup = future_cleanup
        else:
            raise CheckM2DbError("CheckM2 predict CLI does not support a known remove-intermediates option")
    return CheckM2PredictCapabilities(cleanup, _sha256_text(help_text))


def checkm2_database_download_command(executable: Path | str, root: Path | str) -> list[str]:
    return [str(executable), "database", "--download", "--path", str(root), "--no_write_json_db"]


def checkm2_testrun_command(executable: Path | str, database: Path | str, threads: int | str = 1, *, lowmem: bool = False) -> list[str]:
    command = [str(executable)]
    if lowmem:
        command.append("--lowmem")
    command.extend(["testrun", "--threads", str(max(1, int(threads))), "--database_path", str(database)])
    return command


def checkm2_predict_command(executable: Path | str, input_path: Path | str, output_dir: Path | str,
                           database: Path | str, threads: int | str,
                           capabilities: CheckM2PredictCapabilities, *, lowmem: bool = False) -> list[str]:
    command = [str(executable)]
    if lowmem:
        command.append("--lowmem")
    command.extend([
        "predict", "--threads", str(max(1, int(threads))),
        "--input", str(input_path), "--output-directory", str(output_dir),
        "--database_path", str(database), capabilities.cleanup_option, "--force",
    ])
    return command


def checkm2_input_suffix(path: Path) -> str:
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes[-2:] == [".fna", ".gz"]: return ".fna.gz"
    if suffixes[-2:] == [".fa", ".gz"]: return ".fa.gz"
    if suffixes[-2:] == [".fasta", ".gz"]: return ".fasta.gz"
    if suffixes and suffixes[-1] in {".fna", ".fa", ".fasta"}: return suffixes[-1]
    return ".fna.gz" if path.suffix.lower() == ".gz" else ".fna"


def checkm2_named_input_link(assembly: Path, input_dir: Path, isolate: str) -> Path:
    input_dir.mkdir(parents=True, exist_ok=True)
    link = input_dir / f"{safe_name(isolate)}{checkm2_input_suffix(assembly)}"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(assembly.resolve())
    return link


def parse_checkm2_quality_report(path: Path, expected_isolate: str | None = None) -> tuple[float, float]:
    rows = read_tsv(path) if path.is_file() and path.stat().st_size > 0 else []
    if len(rows) != 1:
        raise ValueError(f"CheckM2 report must contain exactly one genome row: {path}")
    row = rows[0]
    observed = (row.get("Name") or row.get("name") or row.get("Bin Id") or row.get("bin_id") or "").strip()
    if expected_isolate and observed and observed != safe_name(expected_isolate):
        raise ValueError(f"CheckM2 report row '{observed}' does not match expected isolate {safe_name(expected_isolate)}")
    completeness = float(row.get("Completeness") or row.get("completeness") or row.get("Completeness_General") or "")
    contamination = float(row.get("Contamination") or row.get("contamination") or "")
    for label, value in (("completeness", completeness), ("contamination", contamination)):
        if not math.isfinite(value) or value < 0 or value > 100:
            raise ValueError(f"CheckM2 {label} is outside 0-100: {value}")
    return completeness, contamination


def companion_python_for_checkm2(executable: Path | str) -> Path:
    exe = Path(executable).expanduser().resolve()
    python = exe.parent / "python"
    return python if python.is_file() else Path("python")


def bundled_test_genome(executable: Path | str) -> Path:
    python = companion_python_for_checkm2(executable)
    script = (
        "import checkm2, pathlib\n"
        "root = pathlib.Path(checkm2.__file__).resolve().parent\n"
        "suffixes = ('.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz')\n"
        "candidates = [p for p in root.rglob('*') if p.is_file() and any(str(p).lower().endswith(s) for s in suffixes)]\n"
        "tests = [p for p in candidates if 'test' in str(p).lower()]\n"
        "pool = tests or candidates\n"
        "print(str(sorted(pool, key=lambda p: (p.stat().st_size, str(p)))[0]) if pool else '')\n"
    )
    try:
        result = subprocess.run([str(python), "-c", script], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        raise CheckM2DbError(f"Could not inspect bundled CheckM2 test genomes with {python}: {error}") from error
    if result.returncode:
        raise CheckM2DbError(f"Could not inspect bundled CheckM2 test genomes: {(result.stderr or result.stdout).strip()}")
    path = Path(result.stdout.strip())
    if not path.is_file():
        raise CheckM2DbError("No bundled CheckM2 test genome could be found in the companion environment")
    return path


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


def checkm2_runtime_signature(cfg: dict[str, str], path: Path | str, executable: Path | str, version: str) -> dict[str, object]:
    database = Path(path).expanduser().resolve()
    resolved_executable = Path(executable).expanduser().resolve()
    env_history = _checkm2_env_history(resolved_executable)
    capabilities = checkm2_predict_capabilities(resolved_executable)
    return {
        "schema": CHECKM2_COMMAND_SCHEMA_VERSION,
        "CHECKM2_DB": str(database),
        "CHECKM2_DB_FILE": _file_signature(database),
        "CHECKM2_EXECUTABLE": str(resolved_executable),
        "CHECKM2_EXECUTABLE_FILE": _file_signature(resolved_executable),
        "CHECKM2_CONDA_HISTORY_FILE": _file_signature(env_history) if env_history else None,
        "CHECKM2_VERSION": version,
        "CHECKM2_LOWMEM": str(truthy(cfg.get("CHECKM2_LOWMEM", "false"))).lower(),
        "CHECKM2_PREDICT_CLEANUP_OPTION": capabilities.cleanup_option,
        "CHECKM2_PREDICT_HELP_SHA256": capabilities.help_sha256,
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
        expected = checkm2_runtime_signature(cfg, path, executable, version)
    except (OSError, ValueError):
        return False
    return recorded.get("status") == "complete" and all(recorded.get(key) == value for key, value in expected.items())


def record_checkm2_runtime_verified(cfg: dict[str, str], path: Path | str,
                                    executable: Path | str, version: str) -> Path:
    marker = checkm2_runtime_marker(cfg)
    atomic_json(marker, {"status": "complete", **checkm2_runtime_signature(cfg, path, executable, version)})
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
        runner(checkm2_database_download_command(executable, root))
        downloaded = find_managed_checkm2_db(root)
        if not downloaded:
            raise CheckM2DbError(f"CheckM2 download completed but {EXPECTED_CHECKM2_DB_NAME} was not found under {root}")
        validate_checkm2_db(downloaded)
        return CheckM2DbResolution(downloaded, "auto_download")
