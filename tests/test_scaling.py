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
from cleangene.util import atomic_json, read_tsv, write_tsv
from cleangene.workers import retained_rows, validate_preflight_input_paths


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


if __name__ == "__main__":
    unittest.main()
