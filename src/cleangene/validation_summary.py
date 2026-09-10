"""Reconciled matrix summaries and decision-set counts, without per-call buffering."""
from __future__ import annotations

import csv
from collections import Counter
from itertools import zip_longest
from pathlib import Path

from .util import write_tsv

SUMMARY_VERSION = 1
COUNT_FIELDS = ["total_calls", "kept", "changed", "changed_0_to_1", "changed_1_to_0",
                "initial_present", "validated_present"]
PCT_FIELDS = ["kept", "changed", "changed_0_to_1", "changed_1_to_0"]
SUMMARY_FIELDS = ["group_id", "isolate_id", *COUNT_FIELDS, *[f"{f}_pct" for f in PCT_FIELDS]]
SUMMARY_FILES = ("gene_call_summary.tsv", "gene_call_summary.per_isolate.tsv",
                 "summary_statistics.txt", "decision_reason_upset.tsv")
PLOT_FILES = ("decision_reason_upset.png", "decision_reason_upset.svg")
COHORT_FILES = ("gene_call_summary.tsv", "gene_call_summary.per_group.tsv",
                "gene_call_summary.per_isolate.tsv", "decision_reason_upset.tsv", *PLOT_FILES)

# Metric membership means input to the deciding rule, not a threshold passed.
STATE_METRICS = {
    "confirmed_present": ("breadth", "identity", "mean_depth"),
    "divergent_variant": ("breadth", "identity", "mean_depth"),
    "possible_truncation": ("breadth", "identity", "mean_depth"),
    "partial_homolog": ("breadth", "identity", "mean_depth"),
    "ambiguous_multimap": ("unique_mapped_reads", "ambiguous_mapped_reads"),
    "not_detected": ("mapped_reads", "breadth"),
    "insufficient_evidence": ("mean_depth", "identity"),
    "confirmed_absent_locus": ("junction_identity", "junction_spanning_alignment", "flank_anchor_length"),
}


def decision_metric_membership(row: dict) -> tuple[tuple[str, ...], str]:
    recorded = str(row.get("decision_metrics") or "")
    if recorded:
        return tuple(sorted(set(recorded.split(";")))), "recorded"
    state = row.get("evidence_state") or row.get("validation_state")
    metrics = STATE_METRICS.get(state, ("unrecorded_evidence",))
    if state == "ambiguous_multimap" and row.get("family_breadth"):
        metrics = (*metrics, "breadth", "family_breadth")
    return tuple(sorted(metrics)), "inferred_from_state" if state in STATE_METRICS else "unavailable"


def _finish_counts(row: dict) -> dict:
    total = row["total_calls"]
    if row["kept"] + row["changed"] != total:
        raise ValueError("Gene-call summary does not reconcile: kept + changed != total")
    if row["changed_0_to_1"] + row["changed_1_to_0"] != row["changed"]:
        raise ValueError("Gene-call summary does not reconcile: directions != changed")
    if row["initial_present"] + row["changed_0_to_1"] - row["changed_1_to_0"] != row["validated_present"]:
        raise ValueError("Gene-call summary does not reconcile: present counts")
    return {**row, **{f"{key}_pct": round(100 * row[key] / total, 6) if total else 0.0 for key in PCT_FIELDS}}


def _total(rows: list[dict], group: str) -> dict:
    return _finish_counts({"group_id": group, "isolate_id": "ALL",
                           **{key: sum(int(row[key]) for row in rows) for key in COUNT_FIELDS}})


def _binary(value: str) -> int:
    if value not in {"0", "1"}:
        raise ValueError(f"Non-binary gene call: {value!r}")
    return int(value)


def write_validation_summary(initial: Path, validated: Path, evidence: Path,
                             outdir: Path, group: str) -> dict:
    """Validate matching matrices/evidence, then publish group and isolate reports.

    Matrix gene order and evidence (gene, isolate) order must agree with reduction.
    Each retained isolate's denominator is every cluster row, including zeros and
    untested calls carried forward. Counts are gene-isolate cells, not unique genes.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    # Remove the success-bearing report first, so a failed rerun cannot look current.
    (outdir / "summary_statistics.txt").unlink(missing_ok=True)
    intersections = Counter()
    prevalence = Counter()
    states = Counter()
    n_genes = absent = 0
    with initial.open(newline="") as a, validated.open(newline="") as b, evidence.open(newline="") as e:
        initial_rows, final_rows, evidence_rows = (csv.DictReader(h, delimiter="\t") for h in (a, b, e))
        if initial_rows.fieldnames != final_rows.fieldnames or not initial_rows.fieldnames or initial_rows.fieldnames[0] != "Gene":
            raise ValueError("Initial and validated matrix headers differ or lack Gene")
        isolates = initial_rows.fieldnames[1:]
        if len(isolates) != len(set(isolates)):
            raise ValueError("Duplicate isolate columns")
        per_isolate = [{"group_id": group, "isolate_id": iso, **dict.fromkeys(COUNT_FIELDS, 0)} for iso in isolates]
        seen = set()
        for old, new in zip_longest(initial_rows, final_rows):
            if old is None or new is None or old["Gene"] != new["Gene"] or old["Gene"] in seen:
                raise ValueError("Matrix gene rows differ, are missing, or contain duplicates")
            gene = old["Gene"]; seen.add(gene); n_genes += 1; present = 0
            for iso, counts in zip(isolates, per_isolate):
                before, after = _binary(old[iso]), _binary(new[iso])
                row = next(evidence_rows, None)
                if row is None or row.get("Gene") != gene or row.get("isolate_id") != iso:
                    raise ValueError(f"Missing or misordered evidence for {gene}/{iso}")
                if _binary(row.get("initial_call")) != before or _binary(row.get("final_call")) != after:
                    raise ValueError(f"Evidence contradicts matrices for {gene}/{iso}")
                counts["total_calls"] += 1
                counts["initial_present"] += before; counts["validated_present"] += after; present += after
                states[row.get("evidence_state") or row.get("validation_state") or "unrecorded"] += 1
                if before == after:
                    counts["kept"] += 1
                else:
                    direction = "0_to_1" if after else "1_to_0"
                    counts["changed"] += 1; counts[f"changed_{direction}"] += 1
                    metrics, provenance = decision_metric_membership(row)
                    reason = row.get("decision_reason") or row.get("evidence_state") or "unrecorded decision"
                    intersections[(reason, metrics, provenance, direction)] += 1
            # Integer comparisons make exact 15%, 95%, and 99% boundaries stable.
            n = len(isolates)
            bucket = "core" if n and present * 100 >= n * 99 else "soft_core" if n and present * 100 >= n * 95 else "shell" if n and present * 100 >= n * 15 else "cloud"
            prevalence[bucket] += 1; absent += int(present == 0)
        if next(evidence_rows, None) is not None:
            raise ValueError("Evidence contains extra gene-isolate calls")
    per_isolate = [_finish_counts(row) for row in per_isolate]
    total = _total(per_isolate, group)
    if total["total_calls"] != n_genes * len(isolates) or sum(prevalence.values()) != n_genes or sum(intersections.values()) != total["changed"]:
        raise ValueError("Summary reconciliation failed")
    write_tsv(outdir / "gene_call_summary.tsv", SUMMARY_FIELDS, [total])
    write_tsv(outdir / "gene_call_summary.per_isolate.tsv", SUMMARY_FIELDS, per_isolate)
    plot_rows = [{"decision_reason": reason, "decision_metrics": ";".join(metrics),
                  "metrics_provenance": provenance, "change": direction.replace("_to_", "->"), "n_calls": count}
                 for (reason, metrics, provenance, direction), count in sorted(intersections.items())]
    write_tsv(outdir / "decision_reason_upset.tsv", ["decision_reason", "decision_metrics", "metrics_provenance", "change", "n_calls"], plot_rows)
    lines = [f"Validated pangenome statistics: {group}",
             "Source: validated_gene_presence_absence.binary.tsv (final calls, including provisional carry-forward)",
             "Validation: PASS — matrices, evidence, totals, directions, and prevalence bins reconcile.",
             f"Isolates\t{len(isolates)}", f"Total genes\t{n_genes}",
             f"Core genes\t(99% <= isolates <= 100%)\t{prevalence['core']}",
             f"Soft core genes\t(95% <= isolates < 99%)\t{prevalence['soft_core']}",
             f"Shell genes\t(15% <= isolates < 95%)\t{prevalence['shell']}",
             f"Cloud genes\t(0% <= isolates < 15%)\t{prevalence['cloud']}",
             f"Absent in all isolates (included in cloud)\t{absent}", "",
             "Gene-call totals (one call per gene cluster per retained isolate):",
             *[f"{key}\t{total[key]}" for key in COUNT_FIELDS],
             *[f"{key}_pct\t{total[f'{key}_pct']:.6f}" for key in PCT_FIELDS], "",
             "Kept means unchanged, including 0->0 and 1->1; it does not mean independently confirmed.",
             "Percentages use all matrix gene clusters per isolate; total percentages use all gene-isolate cells.",
             "Per-isolate counts: gene_call_summary.per_isolate.tsv", "", "Final evidence states:",
             *[f"{state}\t{count}" for state, count in sorted(states.items())], "",
             "UpSet dots identify decision-rule inputs, not passed thresholds or causal attribution.",
             "Legacy metric membership is inferred from evidence state and labeled in the plot data."]
    temporary = outdir / ".summary_statistics.txt.tmp"
    temporary.write_text("\n".join(lines) + "\n"); temporary.replace(outdir / "summary_statistics.txt")
    return total


def combine_validation_summaries(group_dirs: list[Path], outdir: Path) -> None:
    """Cohort totals count calls across groups, without pooling gene identifiers."""
    from .util import read_tsv
    group_rows, isolate_rows = [], []
    intersections = Counter()
    for directory in group_dirs:
        group_rows.extend(read_tsv(directory / "gene_call_summary.tsv"))
        isolate_rows.extend(read_tsv(directory / "gene_call_summary.per_isolate.tsv"))
        for row in read_tsv(directory / "decision_reason_upset.tsv"):
            intersections[(row["decision_reason"], row["decision_metrics"], row["metrics_provenance"], row["change"])] += int(row["n_calls"])
    total = _total(group_rows, "ALL")
    write_tsv(outdir / "gene_call_summary.tsv", SUMMARY_FIELDS, [total])
    write_tsv(outdir / "gene_call_summary.per_group.tsv", SUMMARY_FIELDS, group_rows)
    write_tsv(outdir / "gene_call_summary.per_isolate.tsv", SUMMARY_FIELDS, isolate_rows)
    if sum(intersections.values()) != total["changed"]:
        raise ValueError("Cohort decision counts do not reconcile with changed calls")
    write_tsv(outdir / "decision_reason_upset.tsv", ["decision_reason", "decision_metrics", "metrics_provenance", "change", "n_calls"],
              [[*key, count] for key, count in sorted(intersections.items())])
