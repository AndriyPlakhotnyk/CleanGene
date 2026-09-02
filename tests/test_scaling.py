import contextlib
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from cleangene.cli import main, make_run
from cleangene.defaults import DEFAULTS
from cleangene.qc import prepare_qc_provenance, resolve_threshold_rows
from cleangene.task_store import build_group_task_store, build_isolate_task_store, load_group_task
from cleangene.util import atomic_json, load_json, read_tsv, write_tsv
from cleangene.workers import global_preflight, retained_rows, slurm_controller, validate_preflight_input_paths


class ScalingArchitectureTests(unittest.TestCase):
    def test_launcher_submits_without_checkm2_or_database_preflight(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); manifest=root/"manifest.tsv"
            manifest.write_text("isolate_id\tgroup_id\tR1\tR2\ni1\tg\tmissing_R1.fastq.gz\tmissing_R2.fastq.gz\n")
            with patch("cleangene.cli.resolve_checkm2_executable",side_effect=AssertionError("CheckM2 resolved in launcher")), \
                 patch("cleangene.cli.resolve_checkm2_db",side_effect=AssertionError("CheckM2 DB checked in launcher")), \
                 patch("cleangene.cli.resolve_kraken2_db",create=True,side_effect=AssertionError("Kraken DB checked in launcher")), \
                 patch("cleangene.cli.command_exists",side_effect=AssertionError("tool check ran in launcher")), \
                 patch("cleangene.cli.slurm",return_value="123") as submit, \
                 contextlib.redirect_stdout(StringIO()) as stdout:
                self.assertEqual(main(["run","--manifest",str(manifest),"--analysis-root",str(root),"--run-id","thin"]),0)
            submit.assert_called_once()
            self.assertIn("Controller submitted: 123",stdout.getvalue())
            timings=read_tsv(root/"runs"/"thin"/"logs"/"launcher_timing.tsv")
            self.assertIn("submit_controller",{row["phase"] for row in timings})

    def test_preflight_reports_multiple_missing_inputs(self):
        rows=[
            {"isolate_id":"i1","group_id":"g","R1":"/missing/a_R1.fq.gz","R2":"/missing/a_R2.fq.gz"},
            {"isolate_id":"i2","group_id":"g","R1":"/missing/b_R1.fq.gz","R2":"/missing/b_R2.fq.gz"},
        ]
        with self.assertRaises(SystemExit) as raised:
            validate_preflight_input_paths(rows,{**DEFAULTS,"PREFLIGHT_FILE_CHECK_WORKERS":"2"})
        message=str(raised.exception)
        for expected in ("i1.R1","i1.R2","i2.R1","i2.R2"):
            self.assertIn(expected,message)

    def test_group_store_limits_retained_rows_to_group_members(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            rows=[
                {"isolate_id":"a","group_id":"g1","R1":"r1","R2":"r2"},
                {"isolate_id":"b","group_id":"g2","R1":"r1","R2":"r2"},
                {"isolate_id":"c","group_id":"g1","R1":"r1","R2":"r2"},
            ]
            atomic_json(run/"provenance"/"resolved_config.json",DEFAULTS)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],rows)
            prepare_qc_provenance(run,rows,DEFAULTS)
            build_isolate_task_store(run,rows); build_group_task_store(run,rows)
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g1",2,"small"],["g2",1,"small"]])
            for iso,group in (("a","g1"),("c","g1")):
                sample=run/"results"/"sample_data"/iso
                write_tsv(sample/"qc.tsv",["isolate_id","group_id","excluded","R1","R2","assembly","gff"],[[iso,group,"0","r1","r2","asm","gff"]])
            self.assertEqual(load_group_task(run,0)["isolate_ids"],["a","c"])
            self.assertEqual([row["isolate_id"] for row in retained_rows(run,"g1")],["a","c"])

    def test_qc_profile_resolution_scales_linearly(self):
        rows=[{"isolate_id":f"i{i}","group_id":f"g{i%10}","organism":"Species"} for i in range(1000)]
        with tempfile.TemporaryDirectory() as d:
            profile=Path(d)/"profiles.tsv"
            profile.write_text("scope_type\tscope_value\tqc_min_n50_pass\norganism\tSpecies\t25000\ngroup_id\tg1\t26000\n")
            resolved=resolve_threshold_rows(rows,DEFAULTS,profile)
        self.assertEqual(len(resolved),1000)
        self.assertEqual(resolved[1]["qc_profile_source"],"group_id:g1")

    def test_global_preflight_does_not_run_database_setup(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); run=root/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            r1=root/"r1.fq"; r2=root/"r2.fq"; r1.write_text("@r\nA\n+\nI\n"); r2.write_text("@r\nT\n+\nI\n")
            cfg={**DEFAULTS,"TAXONOMY_MODE":"kraken2","CHECKM2_MODE":"required","CHECKM2_AUTO_DOWNLOAD":"true","PREFLIGHT_FILE_CHECK_WORKERS":"1"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            atomic_json(run/"provenance"/"inputs.json",{"manifest":str(root/"manifest.tsv")})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["i1","g",str(r1),str(r2)]])
            with patch("cleangene.workers.kraken_db_setup",side_effect=AssertionError("Kraken setup ran in preflight")), \
                 patch("cleangene.workers.checkm2_db_setup",side_effect=AssertionError("CheckM2 setup ran in preflight")), \
                 patch("cleangene.workers.run",side_effect=AssertionError("heavy command ran in preflight")):
                payload=global_preflight(run)
            self.assertEqual(payload["status"],"PASS")
            self.assertEqual(payload["Kraken2 setup"],"required")
            self.assertEqual(payload["CheckM2 setup"],"required")
            self.assertFalse((run/"state"/"kraken_db_setup.done.json").exists())
            self.assertFalse((run/"state"/"checkm2_db_setup.done.json").exists())

    def test_controller_routes_database_setup_to_dedicated_resources_before_preprocess(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg={**DEFAULTS,"TAXONOMY_MODE":"kraken2","CHECKM2_MODE":"required","SLURM_POLL_SECONDS":"0","KRAKEN2_DB_CPUS":"8","KRAKEN2_DB_MEM":"64G","KRAKEN2_DB_TIME":"08:00:00","CHECKM2_CPUS":"16","CHECKM2_MEM":"128G","CHECKM2_TIME":"24:00:00","SLURM_CONTROLLER_CPUS":"1","SLURM_CONTROLLER_MEM":"4G","SLURM_CONTROLLER_TIME":"7-00:00:00"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["i1","g","r1","r2"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["g","i1"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g",1,"small"]])
            order=[]
            def fake_single(run_dir,cfg_arg,stage,cpus,mem,time_limit,label):
                order.append((stage,cpus,mem,time_limit))
                atomic_json(run_dir/"state"/f"{stage}.done.json",{"status":"complete"})
                return stage
            with patch("cleangene.workers.global_preflight",side_effect=lambda rd: order.append(("preflight","","","")) or {"status":"PASS"}), \
                 patch("cleangene.workers._run_single_job",side_effect=fake_single), \
                 patch("cleangene.workers.run_resume_maintenance"), \
                 patch("cleangene.workers.user_queue_snapshot",return_value={"total":0,"jobs":{},"entries":[]}), \
                 patch("cleangene.workers.reconcile_preprocess_outputs",return_value={"total":1,"marker_complete":0,"output_recovered":0,"active":0,"incomplete":1,"inconsistent":0}), \
                 patch("cleangene.workers._controller_pipeline",side_effect=lambda *args: order.append(("preprocess","","",""))):
                slurm_controller(run)
            self.assertEqual(order[:4],[
                ("preflight","","",""),
                ("kraken_db_setup","8","64G","08:00:00"),
                ("checkm2_db_setup","16","128G","24:00:00"),
                ("resolve_groups","2","8G","04:00:00"),
            ])
            self.assertEqual(order[4],("preprocess","","",""))

    def test_controller_stops_before_preprocess_if_setup_job_fails(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg={**DEFAULTS,"TAXONOMY_MODE":"off","CHECKM2_MODE":"required","SLURM_POLL_SECONDS":"0","CHECKM2_CPUS":"16","CHECKM2_MEM":"128G","CHECKM2_TIME":"24:00:00"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["i1","g","r1","r2"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["g","i1"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g",1,"small"]])
            with patch("cleangene.workers.global_preflight",return_value={"status":"PASS"}), \
                 patch("cleangene.workers._run_single_job",side_effect=RuntimeError("setup failed")), \
                 patch("cleangene.workers._controller_pipeline") as pipeline:
                with self.assertRaisesRegex(RuntimeError,"setup failed"):
                    slurm_controller(run)
            pipeline.assert_not_called()

    def test_resume_after_incomplete_checkm2_setup_submits_setup_then_pipeline(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg={**DEFAULTS,"TAXONOMY_MODE":"off","CHECKM2_MODE":"required","SLURM_POLL_SECONDS":"0","CHECKM2_CPUS":"16","CHECKM2_MEM":"128G","CHECKM2_TIME":"24:00:00"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["i1","g","r1","r2"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["g","i1"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g",1,"small"]])
            atomic_json(run/"state"/"preprocess"/"i1.done.json",{"status":"complete"})
            order=[]
            with patch("cleangene.workers.global_preflight",side_effect=lambda rd: order.append("preflight") or {"status":"PASS"}), \
                 patch("cleangene.workers._run_single_job",side_effect=lambda *args: order.append(args[2]) or "jid"), \
                 patch("cleangene.workers.run_resume_maintenance"), \
                 patch("cleangene.workers.user_queue_snapshot",return_value={"total":0,"jobs":{},"entries":[]}), \
                 patch("cleangene.workers.reconcile_preprocess_outputs",return_value={"total":1,"marker_complete":1,"output_recovered":0,"active":0,"incomplete":0,"inconsistent":0}), \
                 patch("cleangene.workers._controller_pipeline",side_effect=lambda *args: order.append("pipeline")):
                slurm_controller(run)
            self.assertEqual(order,["preflight","checkm2_db_setup","resolve_groups","pipeline"])

    def test_resume_ignore_checkm2_skips_stale_setup_failure(self):
        with tempfile.TemporaryDirectory() as d:
            from cleangene.cli import main
            run=Path(d)/"runs"/"failed"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir(parents=True)
            cfg={**DEFAULTS,"TAXONOMY_MODE":"off","CHECKM2_MODE":"required","SLURM_POLL_SECONDS":"0"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["i1","g","r1","r2"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["g","i1"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g",1,"small"]])
            atomic_json(run/"state"/"checkm2_db_setup.done.json",{"status":"failed","reason":"old oom"})
            with patch("cleangene.cli.slurm",return_value="123"), \
                 patch("cleangene.cli.resolve_checkm2_executable",side_effect=AssertionError("CheckM2 resolved during ignored resume")), \
                 contextlib.redirect_stdout(StringIO()):
                self.assertEqual(main(["resume","--run-dir",str(run),"--ignore-checkm2"]),0)
            resolved=load_json(run/"provenance"/"resolved_config.json")
            self.assertEqual(resolved["CHECKM2_MODE"],"off")
            self.assertEqual(resolved["CHECKM2_DISABLED_BY_USER"],"true")


if __name__ == "__main__":
    unittest.main()
