from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Iterable

from .qc import THRESHOLD_COLUMNS, THRESHOLD_DEFAULTS, validate_thresholds
from .util import read_tsv

OFFSET_WIDTH = 8
TASK_DIR = Path("state") / "tasks"
ISOLATE_JSONL = TASK_DIR / "isolate_tasks.jsonl"
ISOLATE_INDEX = TASK_DIR / "isolate_tasks.idx"
GROUP_JSONL = TASK_DIR / "group_tasks.jsonl"
GROUP_INDEX = TASK_DIR / "group_tasks.idx"


def isolate_store_paths(run_dir: Path) -> tuple[Path, Path]:
    return run_dir / ISOLATE_JSONL, run_dir / ISOLATE_INDEX


def group_store_paths(run_dir: Path) -> tuple[Path, Path]:
    return run_dir / GROUP_JSONL, run_dir / GROUP_INDEX


def _threshold_map(run_dir: Path) -> dict[str, dict[str, str]]:
    path = run_dir / "provenance" / "qc_thresholds.tsv"
    if not path.is_file():
        return {}
    return {row["isolate_id"]: row for row in read_tsv(path)}


def build_isolate_task_store(run_dir: Path, rows: Iterable[dict[str, str]]) -> int:
    """Write O(1)-addressable isolate task records.

    The JSONL contains complete manifest/task rows plus resolved QC thresholds.
    The index contains unsigned 8-byte offsets into the JSONL, one per task.
    """
    jsonl, index = isolate_store_paths(run_dir)
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    thresholds = _threshold_map(run_dir)
    tmp_jsonl = jsonl.with_suffix(".jsonl.tmp")
    tmp_index = index.with_suffix(".idx.tmp")
    count = 0
    with tmp_jsonl.open("wb") as data, tmp_index.open("wb") as idx:
        for count, row in enumerate(rows, 1):
            threshold_row = {**THRESHOLD_DEFAULTS, **thresholds.get(row["isolate_id"], {})}
            record = {
                **row,
                "task_index": count - 1,
                "qc_profile_source": threshold_row.get("qc_profile_source", "global"),
                "qc_thresholds": {key: threshold_row.get(key, "") for key in THRESHOLD_COLUMNS},
            }
            idx.write(data.tell().to_bytes(OFFSET_WIDTH, "big", signed=False))
            data.write(json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8") + b"\n")
    os.replace(tmp_jsonl, jsonl)
    os.replace(tmp_index, index)
    return count


def _write_indexed_jsonl(jsonl: Path, index: Path, records: Iterable[dict[str, object]]) -> int:
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    tmp_jsonl = jsonl.with_suffix(jsonl.suffix + ".tmp")
    tmp_index = index.with_suffix(index.suffix + ".tmp")
    count = 0
    with tmp_jsonl.open("wb") as data, tmp_index.open("wb") as idx:
        for count, record in enumerate(records, 1):
            idx.write(data.tell().to_bytes(OFFSET_WIDTH, "big", signed=False))
            data.write(json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8") + b"\n")
    os.replace(tmp_jsonl, jsonl)
    os.replace(tmp_index, index)
    return count


def build_group_task_store(run_dir: Path, rows: Iterable[dict[str, str]]) -> int:
    groups: dict[str, dict[str, object]] = {}
    for index, row in enumerate(rows):
        group = row["group_id"]
        record = groups.setdefault(group, {
            "group_id": group,
            "organism": row.get("organism", ""),
            "isolate_indices": [],
            "isolate_ids": [],
            "pangenome_dir": row.get("pangenome_dir", ""),
        })
        record["isolate_indices"].append(index)
        record["isolate_ids"].append(row["isolate_id"])
        if row.get("organism") and not record.get("organism"):
            record["organism"] = row["organism"]
        if row.get("pangenome_dir") and not record.get("pangenome_dir"):
            record["pangenome_dir"] = row["pangenome_dir"]
    jsonl, index = group_store_paths(run_dir)
    return _write_indexed_jsonl(jsonl, index, groups.values())


def task_store_ready(run_dir: Path) -> bool:
    jsonl, index = isolate_store_paths(run_dir)
    return jsonl.is_file() and index.is_file() and index.stat().st_size % OFFSET_WIDTH == 0


def group_store_ready(run_dir: Path) -> bool:
    jsonl, index = group_store_paths(run_dir)
    return jsonl.is_file() and index.is_file() and index.stat().st_size % OFFSET_WIDTH == 0


def _load_indexed_record(jsonl: Path, idx: Path, index: int, label: str) -> dict[str, object]:
    if index < 0:
        raise SystemExit(f"Task index {index} outside {label} task store")
    with idx.open("rb") as handle:
        handle.seek(index * OFFSET_WIDTH)
        raw = handle.read(OFFSET_WIDTH)
    if len(raw) != OFFSET_WIDTH:
        raise SystemExit(f"Task index {index} outside {label} task store")
    offset = int.from_bytes(raw, "big", signed=False)
    with jsonl.open("rb") as handle:
        handle.seek(offset)
        line = handle.readline()
    if not line:
        raise SystemExit(f"Task index {index} has no {label} task record")
    return json.loads(line)


def load_isolate_task(run_dir: Path, index: int) -> dict[str, object]:
    jsonl, idx = isolate_store_paths(run_dir)
    record = _load_indexed_record(jsonl, idx, index, "isolate")
    thresholds = record.get("qc_thresholds", {})
    record["qc_thresholds_resolved"] = validate_thresholds(thresholds, str(record.get("isolate_id", index)))
    return record


def load_group_task(run_dir: Path, index: int) -> dict[str, object]:
    jsonl, idx = group_store_paths(run_dir)
    return _load_indexed_record(jsonl, idx, index, "group")


def migrate_isolate_task_store(run_dir: Path) -> int:
    if task_store_ready(run_dir):
        return 0
    manifest = run_dir / "provenance" / "manifest.tsv"
    tasks = run_dir / "state" / "isolate_tasks.tsv"
    if manifest.is_file():
        rows = read_tsv(manifest)
    elif tasks.is_file():
        rows = read_tsv(tasks)
    else:
        raise SystemExit(f"Cannot build isolate task store; missing {manifest} and {tasks}")
    return build_isolate_task_store(run_dir, rows)


def migrate_group_task_store(run_dir: Path) -> int:
    if group_store_ready(run_dir):
        return 0
    manifest = run_dir / "provenance" / "manifest.tsv"
    if not manifest.is_file():
        raise SystemExit(f"Cannot build group task store; missing {manifest}")
    return build_group_task_store(run_dir, read_tsv(manifest))
