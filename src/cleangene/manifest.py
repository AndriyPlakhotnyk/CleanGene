from __future__ import annotations
import csv
from pathlib import Path
from .util import write_tsv

PANGENOME_COLUMNS = ("pangenome_dir", "panaroo_dir", "pangenome")
BAM_COLUMNS = ("raw_bam", "BAM", "bam")

def _clean(value: str | None) -> str:
    return (value or "").strip()

def load_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"Manifest not found: {path}")
    with path.open(newline="", errors="replace") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise SystemExit(f"Manifest is empty: {path}")
    reader = csv.DictReader(lines, delimiter="\t")
    fields = set(reader.fieldnames or [])
    required = {"isolate_id"}
    missing = sorted(required - fields)
    if missing:
        raise SystemExit("Manifest missing required columns: " + ", ".join(missing))
    rows = []
    seen = set()
    for line, raw in enumerate(reader, 2):
        row = {k: _clean(v) for k, v in raw.items() if k is not None}
        isolate = row.get("isolate_id", "")
        if not isolate:
            raise SystemExit(f"Manifest row {line} has no isolate_id")
        if isolate in seen:
            raise SystemExit(f"Duplicate isolate_id in manifest: {isolate}")
        seen.add(isolate)
        bam = next((row.get(c, "") for c in BAM_COLUMNS if row.get(c, "")), "")
        row["raw_bam"] = bam
        has_fastq = bool(row.get("R1") and row.get("R2"))
        if not has_fastq and not bam:
            raise SystemExit(f"Manifest row {line} must provide R1/R2 or raw_bam")
        if bool(row.get("R1")) != bool(row.get("R2")):
            raise SystemExit(f"Manifest row {line} must provide both R1 and R2")
        if not row.get("group_id"):
            row["group_id"] = row.get("organism") or "default"
        pangenome = next((row.get(c, "") for c in PANGENOME_COLUMNS if row.get(c, "")), "")
        row["pangenome_dir"] = pangenome
        rows.append(row)
    return rows

def groups(rows: list[dict[str, str]]) -> list[str]:
    seen = []
    known = set()
    for row in rows:
        group = row["group_id"]
        if group not in known:
            known.add(group)
            seen.append(group)
    return seen

def write_resolved(path: Path, rows: list[dict[str, str]]) -> None:
    preferred = ["isolate_id", "group_id", "organism", "R1", "R2", "raw_bam", "assembly", "pangenome_dir"]
    extras = sorted({key for row in rows for key in row} - set(preferred))
    fields = [field for field in preferred if any(field in row for row in rows)] + extras
    write_tsv(path, fields, rows)
