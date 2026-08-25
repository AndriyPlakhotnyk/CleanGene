from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cleangene.completion import reconcile_preprocess_outputs, validate_preprocess_completion
from cleangene.defaults import DEFAULTS
from cleangene.qc import QC_OUTPUT_FIELDS
from cleangene.util import atomic_json, read_tsv, write_tsv
from cleangene.cli import main as cli_main
from cleangene.workers import _ActiveBatch, _RollingScheduler, slurm_controller


class CompletionReconciliationTests(unittest.TestCase):
    def make_run(self, root: Path, n: int = 1, cfg: dict[str, str] | None = None) -> tuple[Path, list[dict[str, str]]]:
        run = root / "run"
        (run / "provenance").mkdir(parents=True)
        (run / "state").mkdir()
        merged = {**DEFAULTS, "TAXONOMY_MODE": "off", "PREPROCESS_USE_NODE_LOCAL_SCRATCH": "false", **(cfg or {})}
        atomic_json(run / "provenance" / "resolved_config.json", merged)
        rows = [{"isolate_id": f"BI_{i:04d}", "group_id": "g", "R1": "r1", "R2": "r2"} for i in range(n)]
        write_tsv(run / "provenance" / "manifest.tsv", ["isolate_id", "group_id", "R1", "R2"], rows)
        write_tsv(run / "state" / "isolate_tasks.tsv", ["group_id", "isolate_id"], ([r["group_id"], r["isolate_id"]] for r in rows))
        write_tsv(run / "state" / "group_tasks.tsv", ["group_id", "n_isolates", "group_size_class"], [["g", n, "small"]])
        return run, rows

    def write_qc(self, run: Path, isolate: str, *, group: str = "g", status: str = "PASS", excluded: str = "0",
                 assembly: str | Path | None = None, gff: str | Path | None = None, extra_fields: bool = True,
                 rows: int = 1) -> Path:
        iso_dir = run / "results" / "groups" / group / "01_isolates" / isolate
        iso_dir.mkdir(parents=True, exist_ok=True)
        if assembly is None:
            assembly_path = iso_dir / "assembly" / "contigs.fa"
            assembly_path.parent.mkdir(exist_ok=True)
            assembly_path.write_text(">c\nACGT\n")
            assembly = assembly_path
        if gff is None:
            gff_path = iso_dir / "annotation" / f"{isolate}.gff"
            gff_path.parent.mkdir(exist_ok=True)
            gff_path.write_text("##gff-version 3\n")
            gff = gff_path
        fields = ["isolate_id", "group_id", "excluded", "reason", "R1", "R2", "assembly", "gff"]
        if extra_fields:
            fields += list(QC_OUTPUT_FIELDS)
        row = {
            "isolate_id": isolate, "group_id": group, "excluded": excluded, "reason": "ok",
            "R1": "r1", "R2": "r2", "assembly": str(assembly or ""), "gff": str(gff or ""),
            "PASS/FAIL": status, "Notes": "All evaluated QC criteria passed",
            "trimmed_read_length": "150", "mean_base_quality": "35", "sequencing_coverage": "40",
            "checkm2_completeness": "95", "checkm2_contamination": "1", "qc_profile_source": "global",
        }
        write_tsv(iso_dir / "qc.tsv", fields, [row for _ in range(rows)])
        return iso_dir / "qc.tsv"

    def test_valid_internal_output_recovers_marker(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            qc = self.write_qc(run, "BI_0000")
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            marker = run / "state" / "preprocess" / "BI_0000.done.json"
            self.assertEqual(counts["output_recovered"], 1)
            self.assertTrue(marker.is_file())
            self.assertEqual(read_tsv(run / "state" / "preprocess_reconciliation.tsv")[0]["qc_path"], str(qc))

    def test_terminal_fail_or_excluded_does_not_require_assembly_gff(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000", status="FAIL", excluded="1", assembly="", gff="")
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            self.assertEqual(counts["output_recovered"], 1)

    def test_external_pangenome_and_assembler_off_do_not_require_gff(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            rows[0]["pangenome_dir"] = "/external/panaroo"
            self.write_qc(run, "BI_0000", assembly="", gff="")
            self.assertEqual(reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})["output_recovered"], 1)
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d), cfg={"ASSEMBLER": "off"})
            self.write_qc(run, "BI_0000", assembly="", gff="")
            self.assertEqual(reconcile_preprocess_outputs(run, {**DEFAULTS, "ASSEMBLER": "off"}, rows, {"entries": []})["output_recovered"], 1)

    def test_old_group_and_kraken_pending_output_is_found(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000", group="__kraken_pending__")
            self.assertEqual(reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})["output_recovered"], 1)

    def test_existing_marker_is_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            marker = run / "state" / "preprocess" / "BI_0000.done.json"
            atomic_json(marker, {"status": "custom"})
            self.write_qc(run, "BI_0000")
            reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            self.assertEqual(marker.read_text(), '{\n  "status": "custom"\n}')

    def test_malformed_wrong_legacy_multirow_and_missing_artifacts_are_inconsistent(self):
        cases = [
            {"extra_fields": False},
            {"rows": 2},
            {"assembly": "missing.fa"},
            {"gff": "missing.gff"},
        ]
        for kwargs in cases:
            with self.subTest(kwargs=kwargs), tempfile.TemporaryDirectory() as d:
                run, rows = self.make_run(Path(d))
                self.write_qc(run, "BI_0000", **kwargs)
                with self.assertRaises(SystemExit):
                    reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
                self.assertEqual(read_tsv(run / "state" / "preprocess_reconciliation.tsv")[0]["state"], "inconsistent")
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_9999")
            result = validate_preprocess_completion(run, DEFAULTS, rows[0], [run / "results" / "groups" / "g" / "01_isolates" / "BI_9999" / "qc.tsv"])
            self.assertEqual(result.state, "inconsistent")

    def test_multiple_qc_candidates_are_inconsistent_unless_override(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d), cfg={"RESUME_REPROCESS_INCONSISTENT_PREPROCESS": "true"})
            rows[0]["group_id"] = "new_group"
            self.write_qc(run, "BI_0000", group="old")
            self.write_qc(run, "BI_0000", group="older")
            cfg = {**DEFAULTS, "RESUME_REPROCESS_INCONSISTENT_PREPROCESS": "true"}
            counts = reconcile_preprocess_outputs(run, cfg, rows, {"entries": []})
            self.assertEqual(counts["inconsistent"], 1)
            self.assertFalse((run / "state" / "preprocess" / "BI_0000.done.json").exists())

    def test_expected_group_candidate_wins_over_fallback_duplicates(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000", group="g")
            self.write_qc(run, "BI_0000", group="old")
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            report = read_tsv(run / "state" / "preprocess_reconciliation.tsv")
            self.assertEqual(counts["output_recovered"], 1)
            self.assertIn("/g/01_isolates/", report[0]["qc_path"])

    def test_malformed_marker_is_quarantined_and_recreated_from_valid_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            marker = run / "state" / "preprocess" / "BI_0000.done.json"
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.write_text("{not-json")
            self.write_qc(run, "BI_0000")
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            self.assertEqual(counts["output_recovered"], 1)
            self.assertTrue(marker.is_file())
            self.assertTrue(list(marker.parent.glob("BI_0000.done.json.stale.*")))

    def test_stale_marker_is_quarantined_when_outputs_are_incomplete(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            marker = run / "state" / "preprocess" / "BI_0000.done.json"
            atomic_json(marker, {"status": "complete"})
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            self.assertEqual(counts["incomplete"], 1)
            self.assertFalse(marker.exists())
            self.assertTrue(list(marker.parent.glob("BI_0000.done.json.stale.*")))

    def test_active_old_preprocess_task_is_not_classified_from_partial_qc(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000", extra_fields=False)
            snapshot = {"entries": [{"name": "cg-preprocess", "command": f"python -m cleangene _worker --run-dir {run}", "task_id": "0"}]}
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, snapshot)
            self.assertEqual(counts["active"], 1)
            self.assertFalse((run / "state" / "preprocess" / "BI_0000.done.json").exists())

    def test_acceptance_1000_recovers_900_and_leaves_100_incomplete(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d), n=1000)
            for row in rows[:900]:
                self.write_qc(run, row["isolate_id"])
            counts = reconcile_preprocess_outputs(run, DEFAULTS, rows, {"entries": []})
            self.assertEqual(counts["output_recovered"], 900)
            self.assertEqual(counts["incomplete"], 100)
            report = read_tsv(run / "state" / "preprocess_reconciliation.tsv")
            self.assertEqual(sum(1 for row in report if row["state"] == "incomplete"), 100)

    def test_controller_stops_before_sbatch_on_inconsistent_output(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000", extra_fields=False)
            with patch("cleangene.workers.user_queue_snapshot", return_value={"total": 0, "jobs": {}, "entries": []}), \
                 patch("cleangene.workers.submit_with_qos_retry", side_effect=AssertionError("sbatch called")):
                with self.assertRaises(SystemExit):
                    slurm_controller(run)

    def test_cli_reconcile_preprocess_dry_run_then_apply(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000")
            with patch("cleangene.cli.user_queue_snapshot", return_value={"total": 0, "jobs": {}, "entries": []}):
                self.assertEqual(cli_main(["reconcile-preprocess", "--run-dir", str(run)]), 0)
                self.assertFalse((run / "state" / "preprocess" / "BI_0000.done.json").exists())
                self.assertEqual(cli_main(["reconcile-preprocess", "--run-dir", str(run), "--apply"]), 0)
                self.assertTrue((run / "state" / "preprocess" / "BI_0000.done.json").is_file())

    def test_targeted_post_job_reconciliation_recovers_completed_batch_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d))
            self.write_qc(run, "BI_0000")
            scheduler = _RollingScheduler(run, DEFAULTS)
            scheduler.active["123"] = _ActiveBatch("123", "preprocess", [0], "0", seen=True, missing_polls=1)
            with patch("cleangene.workers.user_queue_snapshot", return_value={"total": 0, "jobs": {}, "entries": []}), \
                 patch("cleangene.workers.assert_jobs_succeeded", return_value=None):
                scheduler.refresh()
            self.assertTrue((run / "state" / "preprocess" / "BI_0000.done.json").is_file())

    def test_cpu_limit_caps_preprocess_submissions_without_affecting_default(self):
        with tempfile.TemporaryDirectory() as d:
            run, rows = self.make_run(Path(d), n=4, cfg={"SLURM_USER_CPU_LIMIT": "16", "SLURM_CPU_HEADROOM": "0"})
            cfg = {**DEFAULTS, "SLURM_USER_CPU_LIMIT": "16", "SLURM_CPU_HEADROOM": "0", "SLURM_ARRAY_CHUNK_SIZE": "10", "SLURM_MAX_OUTSTANDING_CHUNKS": "10", "SLURM_MAX_PARALLEL": "10", "SLURM_USER_JOB_LIMIT": "100", "SLURM_JOB_HEADROOM": "0"}
            scheduler = _RollingScheduler(run, cfg)
            scheduler.snapshot = {"total": 1, "jobs": {}, "entries": [{"cpus": 8}]}
            with patch("cleangene.workers.submit_with_qos_retry", return_value="123"):
                submitted = scheduler.submit_ready("preprocess", list(range(4)), "4", "1G", "1:00:00", "CleanGene preprocess", 100)
            self.assertEqual(submitted, 2)


if __name__ == "__main__":
    unittest.main()
