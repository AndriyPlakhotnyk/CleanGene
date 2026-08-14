from __future__ import annotations
from pathlib import Path
from .defaults import DEFAULTS

def read_env(path: Path | None) -> dict[str, str]:
    values = dict(DEFAULTS)
    if not path: return values
    if not path.is_file(): raise SystemExit(f"Configuration not found: {path}")
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line: continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip('"').strip("'")
    return values

def truthy(value: str | bool | None) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}
