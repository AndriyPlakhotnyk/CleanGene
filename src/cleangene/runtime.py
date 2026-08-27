from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


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
