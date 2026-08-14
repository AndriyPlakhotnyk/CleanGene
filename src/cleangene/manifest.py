from __future__ import annotations
from pathlib import Path
from .util import read_tsv, safe_name, write_tsv

REQUIRED = ("group_id", "isolate_id", "R1", "R2")
OPTIONAL = ("organism", "assembly", "notes")

def load_manifest(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if not rows: raise SystemExit(f"Manifest is empty: {path}")
    fields = set(rows[0])
    missing = [x for x in REQUIRED if x not in fields]
    if missing: raise SystemExit("Manifest missing columns: " + ", ".join(missing))
    seen: set[str] = set(); safe_isolates: dict[str,str] = {}; safe_groups: dict[str,str] = {}
    for i, row in enumerate(rows, 2):
        for col in REQUIRED:
            if not row.get(col, "").strip(): raise SystemExit(f"Manifest row {i} has blank {col}")
        isolate = row["isolate_id"].strip(); group=row["group_id"].strip(); row["isolate_id"]=isolate; row["group_id"]=group
        if isolate in seen: raise SystemExit(f"Duplicate isolate_id: {isolate}")
        seen.add(isolate)
        si=safe_name(isolate); sg=safe_name(group)
        if si in safe_isolates and safe_isolates[si]!=isolate: raise SystemExit(f"Filesystem-safe isolate ID collision: {safe_isolates[si]!r} and {isolate!r} -> {si!r}")
        if sg in safe_groups and safe_groups[sg]!=group: raise SystemExit(f"Filesystem-safe group ID collision: {safe_groups[sg]!r} and {group!r} -> {sg!r}")
        safe_isolates[si]=isolate; safe_groups[sg]=group
        for key in ("R1", "R2"):
            fp = Path(row[key]).expanduser().resolve()
            if not fp.is_file() or fp.stat().st_size == 0: raise SystemExit(f"Missing/non-empty {key} for {isolate}: {fp}")
            row[key]=str(fp)
        assembly = row.get("assembly", "").strip()
        if assembly:
            ap=Path(assembly).expanduser().resolve()
            if not ap.is_file() or ap.stat().st_size==0: raise SystemExit(f"Assembly missing/empty for {isolate}: {ap}")
            row["assembly"]=str(ap)
    return rows

def groups(rows: list[dict[str, str]]) -> list[str]:
    return sorted({row["group_id"] for row in rows})

def write_resolved(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [*REQUIRED, *OPTIONAL]
    write_tsv(path, fields, ([row.get(f, "") for f in fields] for row in rows))
