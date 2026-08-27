from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


class ToolResolutionError(RuntimeError):
    pass


def _validate_executable(path: Path, name: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise ToolResolutionError(f"{name} executable not found: {resolved}")
    if not os.access(resolved, os.X_OK):
        raise ToolResolutionError(f"{name} executable is not executable: {resolved}")
    return resolved


def resolve_executable(
    name: str,
    configured: str = "",
    python_executable: str | Path | None = None,
    companion_environments: tuple[str, ...] = (),
) -> Path:
    if configured.strip():
        try:
            return _validate_executable(Path(configured), name)
        except ToolResolutionError as error:
            raise ToolResolutionError(f"Explicit {name} executable is invalid: {error}")
    found = shutil.which(name)
    if found:
        return _validate_executable(Path(found), name)
    python = Path(python_executable or sys.executable).expanduser().resolve()
    sibling = python.parent / name
    if sibling.exists():
        return _validate_executable(sibling, name)
    for environment in companion_environments:
        companion = python.parent.parent.parent / environment / "bin" / name
        if companion.exists():
            return _validate_executable(companion, name)
    env_name = os.environ.get("CONDA_DEFAULT_ENV", "") or os.environ.get("VIRTUAL_ENV", "") or "<unknown>"
    raise ToolResolutionError(
        f"{name} executable could not be resolved from the CleanGene runtime environment.\n"
        f"  python={python}\n"
        f"  environment={env_name}\n"
        f"Install or update the CleanGene environment so `{name} --version` works, "
        f"or set {name.upper()}_EXECUTABLE to an exact executable path."
    )


def resolve_checkm2_executable(configured: str = "", python_executable: str | Path | None = None) -> Path:
    try:
        return resolve_executable(
            "checkm2",
            configured,
            python_executable,
            companion_environments=("cleangene-checkm2",),
        )
    except ToolResolutionError as error:
        raise ToolResolutionError(
            f"{error}\n"
            "Repair the complete CleanGene installation from its checkout with:\n"
            "  bash scripts/install_or_update.sh --recreate"
        ) from error


def executable_version(executable: Path | str, label: str) -> str:
    command = [str(executable), "--version"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        raise ToolResolutionError(f"{label} executable could not run `--version`: {error}")
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise ToolResolutionError(f"{label} executable failed `--version` with exit status {result.returncode}: {detail}")
    return (result.stdout or result.stderr or "").strip().splitlines()[0] if (result.stdout or result.stderr).strip() else "unknown"
