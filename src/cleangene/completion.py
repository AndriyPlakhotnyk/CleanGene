from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from .config import assembler_mode, truthy
from .qc import QC_OUTPUT_FIELDS
from .util import read_tsv, safe_name, touch_done, write_tsv

CompletionState = Literal["complete", "incomplete", "inconsistent", "active", "marker_complete", "output_recovered"]

CORE_QC_FIELDS = ("isolate_id", "group_id", "excluded", "reason", "R1", "R2", "assembly", "gff")
VALID_QC_STATUS = {"PASS", "WARNING", "FAIL"}


@dataclass
class CompletionResult:
    state: CompletionState
    reason: str
    isolate_id: str = ""
    group_id: str = ""
    qc_path: Path | None = None
    marker_path: Path | None = None
    qc_status: str = ""
    excluded: bool = False
    assembly_path: str = ""
    gff_path: str = ""
    marker_payload: dict[str, object] | None = None


def marker_path(run_dir: Path, isolate_id: str) -> Path:
    return run_dir / "state" / "preprocess" / f"{safe_name(isolate_id)}.done.json"


def scan_preprocess_qc_outputs(run_dir: Path) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    for qc in (run_dir / "results" / "sample_data").glob("*/qc.tsv"):
        isolate = qc.parent.name
        index.setdefault(isolate, []).append(qc)
    for qc in (run_dir / "results" / "groups").glob("*/01_isolates/*/qc.tsv"):
        isolate = qc.parent.name
        index.setdefault(isolate, []).append(qc)
    return index


def find_isolate_qc_candidates(run_dir: Path, isolate_id: str, group_id: str = "") -> list[Path]:
    safe = safe_name(isolate_id)
    sample = run_dir / "results" / "sample_data" / safe / "qc.tsv"
    if sample.is_file():
        return [sample]
    if group_id:
        expected = run_dir / "results" / "groups" / safe_name(group_id) / "01_isolates" / safe / "qc.tsv"
        if expected.is_file():
            return [expected]
    return sorted((run_dir / "results" / "groups").glob(f"*/01_isolates/{safe}/qc.tsv"))


def _indexed_candidates(run_dir: Path, task: dict[str, str], qc_index: dict[str, list[Path]]) -> list[Path]:
    safe = safe_name(task["isolate_id"])
    group = task.get("group_id", "")
    sample = run_dir / "results" / "sample_data" / safe / "qc.tsv"
    if sample.is_file():
        return [sample]
    if group:
        expected = run_dir / "results" / "groups" / safe_name(group) / "01_isolates" / safe / "qc.tsv"
        if expected.is_file():
            return [expected]
    return sorted(qc_index.get(safe, []))


def _resolve_artifact(run_dir: Path, qc_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    run_relative = run_dir / path
    if run_relative.exists():
        return run_relative
    return qc_path.parent / path


def _read_single_qc(qc_path: Path) -> tuple[dict[str, str] | None, str]:
    try:
        rows = read_tsv(qc_path)
    except (OSError, UnicodeError, ValueError) as error:
        return None, f"qc_unreadable:{error}"
    if not rows:
        return None, "qc_empty"
    if len(rows) != 1:
        return None, f"qc_row_count={len(rows)}"
    return rows[0], ""


def validate_preprocess_completion(
    run_dir: Path,
    cfg: dict[str, str],
    task: dict[str, str],
    candidates: list[Path],
) -> CompletionResult:
    isolate = task["isolate_id"]
    if not candidates:
        return CompletionResult("incomplete", "missing_marker_and_qc", isolate_id=isolate, group_id=task.get("group_id", ""), marker_path=marker_path(run_dir, isolate))
    if len(candidates) > 1:
        return CompletionResult("inconsistent", "multiple_qc_candidates", isolate, task.get("group_id", ""), candidates[0], marker_path(run_dir, isolate))
    qc_path = candidates[0]
    row, error = _read_single_qc(qc_path)
    if row is None:
        return CompletionResult("inconsistent", error, isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate))
    header = set(row)
    missing = sorted(set(CORE_QC_FIELDS).union(QC_OUTPUT_FIELDS) - header)
    if missing:
        return CompletionResult("inconsistent", "qc_missing_fields:" + ",".join(missing), isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate))
    if row.get("isolate_id") != isolate:
        return CompletionResult("inconsistent", f"qc_wrong_isolate:{row.get('isolate_id','')}", isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate))
    qc_status = row.get("PASS/FAIL", "")
    if qc_status not in VALID_QC_STATUS:
        return CompletionResult("inconsistent", f"invalid_qc_status:{qc_status}", isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate), qc_status)

    excluded = truthy(row.get("excluded", "false")) or qc_status == "FAIL"
    if excluded and not row.get("reason", "").strip():
        return CompletionResult("inconsistent", "excluded_missing_reason", isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate), qc_status, excluded)
    external = bool(task.get("pangenome_dir", "").strip())
    assembly_mode_off = assembler_mode(cfg) == "off"
    assembly_path = row.get("assembly", "").strip()
    gff_path = row.get("gff", "").strip()
    if not excluded and not external and not assembly_mode_off:
        for field in ("assembly", "gff"):
            value = row.get(field, "").strip()
            if not value:
                return CompletionResult("inconsistent", f"missing_{field}", isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate), qc_status, excluded, assembly_path, gff_path)
            artifact = _resolve_artifact(run_dir, qc_path, value)
            if not artifact.is_file() or artifact.stat().st_size == 0:
                return CompletionResult("inconsistent", f"{field}_not_readable:{value}", isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate), qc_status, excluded, assembly_path, gff_path)

    payload = {
        "status": "complete",
        "reconciled_from_outputs": True,
        "isolate_id": isolate,
        "group_id": task.get("group_id", ""),
        "qc_status": qc_status,
        "excluded": excluded,
        "reason": row.get("reason", ""),
        "qc_tsv": str(qc_path),
        "assembly": assembly_path,
        "gff": row.get("gff", "").strip(),
        "reconciled_at": datetime.now(timezone.utc).isoformat(),
    }
    if external:
        payload["external_pangenome"] = task.get("pangenome_dir", "").strip()
    return CompletionResult("complete", "valid_terminal_qc", isolate, task.get("group_id", ""), qc_path, marker_path(run_dir, isolate), qc_status, excluded, assembly_path, gff_path, payload)


def _read_marker(marker: Path) -> tuple[dict[str, object] | None, str]:
    try:
        payload = json.loads(marker.read_text())
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return None, f"marker_malformed:{error}"
    if not isinstance(payload, dict) or payload.get("status") not in {"complete", "custom"}:
        return None, "marker_malformed:invalid_payload"
    return payload, ""


def _quarantine_marker(marker: Path, reason: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target = marker.with_name(f"{marker.name}.stale.{stamp}")
    counter = 1
    while target.exists():
        target = marker.with_name(f"{marker.name}.stale.{stamp}.{counter}")
        counter += 1
    os.replace(marker, target)
    return target


def active_preprocess_indices(run_dir: Path, snapshot: dict[str, object]) -> set[int]:
    active: set[int] = set()
    run_arg = str(run_dir)
    for entry in snapshot.get("entries", []):
        name = str(entry.get("name", ""))
        command = str(entry.get("command", ""))
        if name != "cg-preprocess" or run_arg not in command:
            continue
        try:
            active.add(int(str(entry.get("task_id", ""))))
        except ValueError:
            continue
    return active


def reconcile_preprocess_outputs(
    run_dir: Path,
    cfg: dict[str, str],
    isolate_rows: list[dict[str, str]],
    snapshot: dict[str, object] | None = None,
    *,
    indices: list[int] | None = None,
    apply: bool = True,
) -> dict[str, int]:
    qc_index = scan_preprocess_qc_outputs(run_dir)
    active = active_preprocess_indices(run_dir, snapshot or {})
    rows = []
    counts = {key: 0 for key in ("total", "marker_complete", "output_recovered", "active", "incomplete", "inconsistent")}
    selected = [(i, isolate_rows[i]) for i in (indices if indices is not None else range(len(isolate_rows))) if 0 <= i < len(isolate_rows)]
    counts["total"] = len(selected)
    inconsistent: list[str] = []
    report = run_dir / "state" / "preprocess_reconciliation.tsv"

    for index, task in selected:
        isolate = task["isolate_id"]
        marker = marker_path(run_dir, isolate)
        candidates = _indexed_candidates(run_dir, task, qc_index)
        result = validate_preprocess_completion(run_dir, cfg, task, candidates)
        qc_path = str(result.qc_path or "")
        marker_ok = False
        marker_error = ""
        if marker.is_file():
            marker_payload, marker_error = _read_marker(marker)
            marker_ok = marker_payload is not None
        if marker.is_file() and marker_ok and result.state == "complete":
            state = "marker_complete"; reason = "existing_marker"; counts["marker_complete"] += 1
        elif index in active:
            state = "active"; reason = "live_preprocess_task"; counts["active"] += 1
        else:
            if result.state == "complete":
                if marker.is_file() and not marker_ok:
                    if apply:
                        _quarantine_marker(marker, marker_error)
                    reason = f"{marker_error};{result.reason}"
                else:
                    reason = result.reason
                if apply:
                    touch_done(marker, result.marker_payload or {})
                state = "output_recovered"; counts["output_recovered"] += 1
            elif result.state == "incomplete":
                if marker.is_file() and apply:
                    _quarantine_marker(marker, "stale_marker_incomplete_outputs")
                state = "incomplete"; reason = ("stale_marker_incomplete_outputs;" if marker.is_file() else "") + result.reason; counts["incomplete"] += 1
            else:
                if marker.is_file() and apply:
                    _quarantine_marker(marker, "stale_marker_inconsistent_outputs")
                state = "inconsistent"; reason = ("stale_marker_inconsistent_outputs;" if marker.is_file() else "") + result.reason; counts["inconsistent"] += 1; inconsistent.append(isolate)
        will_submit = state == "incomplete" or (state == "inconsistent" and truthy(cfg.get("RESUME_REPROCESS_INCONSISTENT_PREPROCESS", "false")))
        rows.append({"index": index, "isolate_id": isolate, "group_id": task.get("group_id", ""), "state": state,
            "reason": reason, "qc_path": qc_path, "marker_path": str(marker), "will_submit": "1" if will_submit else "0"})

    write_tsv(report, ["index", "isolate_id", "group_id", "state", "reason", "qc_path", "marker_path", "will_submit"], rows)
    if inconsistent and not truthy(cfg.get("RESUME_REPROCESS_INCONSISTENT_PREPROCESS", "false")):
        preview = ", ".join(inconsistent[:10])
        raise SystemExit(f"Inconsistent preprocess outputs for {len(inconsistent)} isolate(s): {preview}. See {report}. Set RESUME_REPROCESS_INCONSISTENT_PREPROCESS=true to reprocess them without deleting old outputs.")
    return counts
