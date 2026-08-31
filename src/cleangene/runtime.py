from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from .util import atomic_json, load_json, sha256


def source_checkout_root(path: Path | str | None) -> Path | None:
    if path is None:
        return None
    current = Path(path).expanduser().resolve()
    if current.is_file():
        current = current.parent
    for parent in (current, *current.parents):
        if (parent / "pyproject.toml").is_file() and (parent / "src" / "cleangene").is_dir():
            return parent
    return None


def cleangene_module_path() -> Path:
    import cleangene
    return Path(cleangene.__file__).resolve()


def cleangene_project_root() -> Path | None:
    return source_checkout_root(cleangene_module_path())


def cleangene_commit(root: Path | None = None) -> str:
    root = root or cleangene_project_root()
    if root is None:
        return "unknown"
    try:
        result = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, check=True)
        return result.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"

def cleangene_dirty(root: Path | None = None) -> str:
    root = root or cleangene_project_root()
    if root is None:
        return "unknown"
    try:
        result = subprocess.run(["git", "status", "--porcelain"], cwd=root, capture_output=True, text=True, check=True)
        return "dirty" if result.stdout.strip() else "clean"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _configured_root_outside_checkouts(cfg: dict[str, str], checkouts: list[Path]) -> bool:
    for key in ("CLEANGENE_DATABASE_ROOT", "KRAKEN2_DATABASE_ROOT", "CHECKM2_DATABASE_ROOT"):
        value = cfg.get(key, "").strip()
        if not value:
            continue
        root = Path(value).expanduser().resolve()
        if all(root != checkout and checkout not in root.parents for checkout in checkouts):
            return True
    return False


def assert_config_matches_runtime(config: Path | None, cfg: dict[str, str]) -> None:
    if not config:
        return
    active = cleangene_project_root()
    config_root = source_checkout_root(config)
    if not active or not config_root or active == config_root:
        return
    if _configured_root_outside_checkouts(cfg, [active, config_root]):
        return
    raise SystemExit(
        "CleanGene checkout mismatch:\n"
        f"  executable is loading CleanGene from: {active}\n"
        f"  --config belongs to: {config_root}\n"
        "Install or update the intended checkout with:\n"
        f"  cd {config_root}\n"
        "  bash scripts/install_or_update.sh\n"
        "  conda activate cleangene"
    )


@dataclass(frozen=True)
class RuntimeIdentity:
    python: Path
    module_path: Path
    project_root: Path | None
    commit: str
    config: Path | None


def runtime_identity(config: Path | None = None) -> RuntimeIdentity:
    root = cleangene_project_root()
    return RuntimeIdentity(
        python=Path(sys.executable).resolve(),
        module_path=cleangene_module_path(),
        project_root=root,
        commit=cleangene_commit(root),
        config=config.expanduser().resolve() if config else None,
    )


def print_runtime_identity(config: Path | None = None) -> None:
    identity = runtime_identity(config)
    print("CleanGene runtime:", flush=True)
    print(f"  version/commit={identity.commit}", flush=True)
    print(f"  project_root={identity.project_root or '<installed>'}", flush=True)
    print(f"  python={identity.python}", flush=True)
    print(f"  config={identity.config or '<defaults>'}", flush=True)


def runtime_provenance(cfg: dict[str, str] | None = None) -> dict[str, object]:
    cfg = cfg or {}
    root = cleangene_project_root()
    module = cleangene_module_path()
    workers = module.parent / "workers.py"
    checkm2 = module.parent / "checkm2.py"
    payload = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": cleangene_commit(root),
        "git_dirty_status": cleangene_dirty(root),
        "cleangene_package_path": str(module.parent),
        "python_executable": str(Path(sys.executable).resolve()),
        "project_root": str(root) if root else "",
        "workers_py_sha256": sha256(workers) if workers.is_file() else "",
        "checkm2_py_sha256": sha256(checkm2) if checkm2.is_file() else "",
        "checkm2_executable": cfg.get("CHECKM2_EXECUTABLE", ""),
        "checkm2_version": cfg.get("CHECKM2_VERSION", ""),
        "checkm2_database": cfg.get("CHECKM2_DB", ""),
    }
    try:
        from .checkm2 import CHECKM2_COMMAND_SCHEMA_VERSION
        payload["checkm2_command_schema_version"] = CHECKM2_COMMAND_SCHEMA_VERSION
    except Exception:
        payload["checkm2_command_schema_version"] = "unknown"
    return payload


def record_runtime_provenance(run_dir: Path, cfg: dict[str, str]) -> Path:
    path = run_dir / "provenance" / "runtime.json"
    if path.is_file():
        history = run_dir / "provenance" / f"runtime.pre_resume.{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        path.replace(history)
        legacy = run_dir / "provenance" / "runtime.pre_resume.json"
        if not legacy.is_file():
            atomic_json(legacy, load_json(history))
    atomic_json(path, runtime_provenance(cfg))
    return path


def verify_worker_runtime(run_dir: Path) -> None:
    expected_path = run_dir / "provenance" / "runtime.json"
    if not expected_path.is_file():
        return
    try:
        expected = load_json(expected_path)
    except (OSError, ValueError):
        return
    current = runtime_provenance({})
    mismatches = []
    for key in ("git_commit", "cleangene_package_path", "python_executable", "workers_py_sha256", "checkm2_py_sha256", "checkm2_command_schema_version"):
        if str(expected.get(key, "")) != str(current.get(key, "")):
            mismatches.append(key)
    if mismatches:
        raise SystemExit(
            "CleanGene runtime mismatch:\n"
            f"  controller source commit = {expected.get('git_commit','unknown')}\n"
            f"  worker source commit = {current.get('git_commit','unknown')}\n"
            f"  controller package path = {expected.get('cleangene_package_path','unknown')}\n"
            f"  worker package path = {current.get('cleangene_package_path','unknown')}\n"
            f"  mismatched fields = {', '.join(mismatches)}\n"
            "Update/reinstall CleanGene and resume the run."
        )
