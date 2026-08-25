from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .config import assembler_mode, truthy
from .qc import QC_OUTPUT_FIELDS
from .util import read_tsv, safe_name, touch_done, write_tsv

CompletionState = Literal["complete", "incomplete", "inconsistent", "active", "marker_complete"]

CORE_QC_FIELDS = ("isolate_id", "group_id", "excluded", "reason", "R1", "R2", "assembly", "gff")
VALID_QC_STATUS = {"PASS", "WARNING", "FAIL"}


@dataclass
class CompletionResult:
    state: CompletionState
    reason: str
    qc_path: Path | None = None
    marker_payload: dict[str, object] | None = None


def marker_path(run_dir: Path, isolate_id: str) -> Path:
    return run_dir / "state" / "preprocess" / f"{safe_name(isolate_id)}.done.json"


def scan_preprocess_qc_outputs(run_dir: Path) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    for qc in (run_dir / "results" / "groups").glob("*/01_isolates/*/qc.tsv"):
        isolate = qc.parent.name
        index.setdefault(isolate, []).append(qc)
    return index


def find_isolate_qc_candidates(run_dir: Path, isolate_id: str) -> list[Path]:
    safe = safe_name(isolate_id)
    return list((run_dir / "results" / "groups").glob(f"*/01_isolates/{safe}/qc.tsv"))


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
        return CompletionResult("incomplete", "missing_marker_and_qc")
    if len(candidates) > 1:
        return CompletionResult("inconsistent", "multiple_qc_candidates", candidates[0])
    qc_path = candidates[0]
    row, error = _read_single_qc(qc_path)
    if row is None:
        return CompletionResult("inconsistent", error, qc_path)
    header = set(row)
    missing = sorted(set(CORE_QC_FIELDS).union(QC_OUTPUT_FIELDS) - header)
    if missing:
        return CompletionResult("inconsistent", "qc_missing_fields:" + ",".join(missing), qc_path)
    if row.get("isolate_id") != isolate:
        return CompletionResult("inconsistent", f"qc_wrong_isolate:{row.get('isolate_id','')}", qc_path)
    qc_status = row.get("PASS/FAIL", "")
    if qc_status not in VALID_QC_STATUS:
        return CompletionResult("inconsistent", f"invalid_qc_status:{qc_status}", qc_path)

    excluded = truthy(row.get("excluded", "false")) or qc_status == "FAIL"
    external = bool(task.get("pangenome_dir", "").strip())
    assembly_mode_off = assembler_mode(cfg) == "off"
    if not excluded and not external and not assembly_mode_off:
        for field in ("assembly", "gff"):
            value = row.get(field, "").strip()
            if not value:
                return CompletionResult("inconsistent", f"missing_{field}", qc_path)
            artifact = _resolve_artifact(run_dir, qc_path, value)
            if not artifact.is_file() or artifact.stat().st_size == 0:
                return CompletionResult("inconsistent", f"{field}_not_readable:{value}", qc_path)

    payload = {
        "status": "complete",
        "reconciled_from_outputs": True,
        "isolate_id": isolate,
        "qc_status": qc_status,
        "excluded": excluded,
        "reason": row.get("reason", ""),
        "qc_tsv": str(qc_path),
        "gff": row.get("gff", "").strip(),
    }
    if external:
        payload["external_pangenome"] = task.get("pangenome_dir", "").strip()
    return CompletionResult("complete", "valid_terminal_qc", qc_path, payload)


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
) -> dict[str, int]:
    qc_index = scan_preprocess_qc_outputs(run_dir)
    active = active_preprocess_indices(run_dir, snapshot or {})
    rows = []
    counts = {key: 0 for key in ("total", "marker_complete", "output_recovered", "active", "incomplete", "inconsistent")}
    counts["total"] = len(isolate_rows)
    inconsistent: list[str] = []
    report = run_dir / "state" / "preprocess" / "reconciliation.tsv"

    for index, task in enumerate(isolate_rows):
        isolate = task["isolate_id"]
        marker = marker_path(run_dir, isolate)
        qc_path = ""
        if marker.is_file():
            state = "marker_complete"; reason = "existing_marker"; counts["marker_complete"] += 1
        elif index in active:
            state = "active"; reason = "live_preprocess_task"; counts["active"] += 1
        else:
            result = validate_preprocess_completion(run_dir, cfg, task, qc_index.get(safe_name(isolate), []))
            qc_path = str(result.qc_path or "")
            if result.state == "complete":
                touch_done(marker, result.marker_payload or {})
                state = "output_recovered"; reason = result.reason; counts["output_recovered"] += 1
            elif result.state == "incomplete":
                state = "incomplete"; reason = result.reason; counts["incomplete"] += 1
            else:
                state = "inconsistent"; reason = result.reason; counts["inconsistent"] += 1; inconsistent.append(isolate)
        rows.append({"index": index, "isolate_id": isolate, "group_id": task.get("group_id", ""), "state": state,
            "reason": reason, "qc_path": qc_path, "marker_path": str(marker)})

    write_tsv(report, ["index", "isolate_id", "group_id", "state", "reason", "qc_path", "marker_path"], rows)
    if inconsistent and not truthy(cfg.get("RESUME_REPROCESS_INCONSISTENT_PREPROCESS", "false")):
        preview = ", ".join(inconsistent[:10])
        raise SystemExit(f"Inconsistent preprocess outputs for {len(inconsistent)} isolate(s): {preview}. See {report}. Set RESUME_REPROCESS_INCONSISTENT_PREPROCESS=true to reprocess them without deleting old outputs.")
    return counts
