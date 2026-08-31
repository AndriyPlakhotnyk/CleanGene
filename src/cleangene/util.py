from __future__ import annotations
import csv, hashlib, json, os, re, shutil, subprocess, tempfile
from pathlib import Path
from typing import Iterable, Sequence

SAFE_RE = re.compile(r"[^A-Za-z0-9_.=-]+")

def safe_name(value: str) -> str:
    return SAFE_RE.sub("_", value.strip()) or "unnamed"

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Sequence[object] | dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        first = None
        rows = iter(rows)
        try: first = next(rows)
        except StopIteration: first = None
        if isinstance(first, dict):
            writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            if first is not None: writer.writerow(first)
            writer.writerows(rows)
        else:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            if first is not None: writer.writerow(first)
            writer.writerows(rows)

def atomic_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle: json.dump(data, handle, indent=2, sort_keys=True)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)

def load_json(path: Path):
    with path.open() as handle: return json.load(handle)

def run(cmd: Sequence[str], *, stdout: Path | None = None, stderr: Path | None = None,
        cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    out = stdout.open("w") if stdout else None
    err = stderr.open("w") if stderr else None
    try: subprocess.run(list(cmd), check=True, stdout=out, stderr=err, cwd=cwd, env=env)
    finally:
        if out: out.close()
        if err: err.close()

def command_exists(name: str) -> bool:
    return shutil.which(name) is not None

def sha256(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block): h.update(chunk)
    return h.hexdigest()

def touch_done(path: Path, payload: dict[str, object] | None = None) -> None:
    atomic_json(path, payload or {"status": "complete"})
