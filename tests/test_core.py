import tempfile, unittest
from pathlib import Path
from io import StringIO
import contextlib
from cleangene.manifest import groups, load_manifest
from cleangene.fasta import assembly_metrics
from cleangene.pangenome import normalize_panaroo, present, select_rows
from cleangene.slurm import array_task_count, available_slots, submit_with_qos_retry, user_job_count, array_chunks, sbatch_cmd, submit, submit_chunked_arrays
from cleangene.util import atomic_json, write_tsv
from cleangene.workers import _run_array_stage, _wait_jobs, controller_downstream, ensure_kraken2_db, manifest_pangenome_dir, manifest_row_for_task, panaroo, parse_kraken_report, slurm_controller, task_row
from cleangene.cli import make_run, slurm
from unittest.mock import patch
import subprocess

class CleanGeneCoreTests(unittest.TestCase):
    def test_present(self):
        self.assertEqual(present(""),0); self.assertEqual(present("0"),0); self.assertEqual(present("abc_1"),1)

    def test_validation_selection(self):
        isolates=["a","b","c"]
        rows=[{"Gene":"g1","a":1,"b":1,"c":1},{"Gene":"g2","a":1,"b":0,"c":0},{"Gene":"g3","a":0,"b":0,"c":0}]
        self.assertEqual([r["Gene"] for r in select_rows(rows,isolates,"differential",0.95)],["g2"])
        self.assertEqual([r["Gene"] for r in select_rows(rows,isolates,"all",0.95)],["g1","g2","g3"])

    def test_assembly_metrics(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"x.fa"; p.write_text(">a\nAAAA\n>b\nGGNN\n")
            m=assembly_metrics(p); self.assertEqual(m["assembly_length"],8); self.assertEqual(m["contigs"],2); self.assertEqual(m["ambiguous_bases"],2)

    def test_kraken_contamination(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"k.tsv"; p.write_text("70.0\t700\t700\tS\t1\t  Expected species\n6.0\t60\t60\tS\t2\t  Other species\n")
            top, contamination, _=parse_kraken_report(p,"Expected species")
            self.assertEqual(top,"Expected species"); self.assertAlmostEqual(contamination,6.0)

    def test_manifest_infers_group_and_pangenome(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"manifest.tsv"
            p.write_text("# comment\nisolate_id\torganism\tR1\tR2\tpangenome\niso1\tSpecies A\tr1.fq\tr2.fq\t/pan\n")
            rows=load_manifest(p)
            self.assertEqual(rows[0]["group_id"],"Species A")
            self.assertEqual(rows[0]["pangenome_dir"],"/pan")
            self.assertEqual(groups(rows),["Species A"])

    def test_manifest_accepts_raw_bam_without_fastq(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"manifest.tsv"
            p.write_text("isolate_id\tgroup_id\traw_bam\niso1\tg1\treads.bam\n")
            rows=load_manifest(p)
            self.assertEqual(rows[0]["raw_bam"],"reads.bam")
            self.assertEqual(rows[0]["group_id"],"g1")

    def test_manifest_marks_missing_organism_for_kraken_grouping(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"manifest.tsv"
            p.write_text("isolate_id\tR1\tR2\niso1\tr1.fq\tr2.fq\n")
            rows=load_manifest(p)
            self.assertEqual(rows[0]["group_id"],"__kraken_pending__")
            self.assertEqual(rows[0]["grouping_source"],"kraken_pending")

    def test_manifest_tolerates_bom_spaces_and_csv(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"manifest.csv"
            p.write_text("\ufeff isolate_id , organism , R1 , R2 \niso1,Species A,r1.fq,r2.fq\n")
            rows=load_manifest(p)
            self.assertEqual(rows[0]["isolate_id"],"iso1")
            self.assertEqual(rows[0]["group_id"],"Species A")

    def test_worker_task_preserves_fastq_paths_from_four_column_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fq"; r2=root/"r2.fq"
            r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n")
            manifest=root/"manifest.tsv"
            manifest.write_text(f"isolate_id\torganism\tR1\tR2\niso1\tSpecies A\t{r1}\t{r2}\n")
            run=make_run(manifest,root,{"SLURM_ACCOUNT":"","SLURM_PARTITION":""},"r")
            rows=load_manifest(run/"provenance"/"manifest.tsv")
            row=manifest_row_for_task(task_row(run,"isolate",0),rows)
            self.assertEqual(row["isolate_id"],"iso1")
            self.assertEqual(row["organism"],"Species A")
            self.assertEqual(row["R1"],str(r1))
            self.assertEqual(row["R2"],str(r2))

    def test_normalize_panaroo_accepts_safe_isolate_names(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"gene_presence_absence.csv"
            p.write_text("Gene,iso_1,iso_2\nabc,locus1,\ndef,,locus2\n")
            rows=normalize_panaroo(p,["iso 1","iso 2"])
            self.assertEqual(rows,[{"Gene":"abc","iso 1":1,"iso 2":0},{"Gene":"def","iso 1":0,"iso 2":1}])

    def test_manifest_pangenome_dir_requires_one_valid_panaroo_dir(self):
        with tempfile.TemporaryDirectory() as d:
            pan=Path(d)/"panaroo"
            pan.mkdir()
            (pan/"gene_presence_absence.csv").write_text("Gene,iso1,iso2\nabc,a,b\n")
            rows=[{"group_id":"g","pangenome_dir":str(pan)},{"group_id":"g","pangenome_dir":str(pan)}]
            self.assertEqual(manifest_pangenome_dir(rows,"g"),pan.resolve())

    def test_sbatch_array_uses_configured_parallel_value(self):
        cmd=sbatch_cmd(name="x",wrap="echo ok",cpus="1",mem="1G",time="00:05:00",array=array_chunks(range(28294),500,"100")[0])
        self.assertIn("--array",cmd)
        self.assertEqual(cmd[cmd.index("--array")+1],"0-499%100")

    def test_chunked_arrays_chain_dependencies(self):
        seen=[]
        def build(array, dep):
            seen.append((array,dep)); return ["sbatch","--array",array]
        with contextlib.redirect_stdout(StringIO()):
            ids=submit_chunked_arrays(build,array_chunks(range(1200),500,"100"),True,1,"dbjob")
        self.assertEqual(ids,["DRYRUN","DRYRUN","DRYRUN"])
        self.assertEqual(seen,[("0-499%100","dbjob"),("500-999%100","DRYRUN"),("1000-1199%100","DRYRUN")])

    def test_submit_reports_sbatch_stderr(self):
        failed=subprocess.CompletedProcess(["sbatch"],1,stdout="",stderr="Batch job submission failed: throttled")
        with patch("cleangene.slurm.subprocess.run",return_value=failed):
            with self.assertRaisesRegex(RuntimeError,"throttled"):
                submit(["sbatch","--array","0-10%100"],False)

    def test_slurm_dry_run_uses_kraken_db_dependency_and_chunks(self):
        with tempfile.TemporaryDirectory() as d:
            manifest=Path(d)/"m.tsv"
            r1=Path(d)/"r1.fq"; r2=Path(d)/"r2.fq"
            r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n")
            manifest.write_text("isolate_id\torganism\tR1\tR2\n" + "\n".join(f"i{x}\tSpecies\t{r1}\t{r2}" for x in range(3)) + "\n")
            cfg={"SLURM_ACCOUNT":"","SLURM_PARTITION":"","SLURM_CONTROLLER_CPUS":"1","SLURM_CONTROLLER_MEM":"1G","SLURM_CONTROLLER_TIME":"1:00:00","TAXONOMY_MODE":"kraken2"}
            run=make_run(manifest,Path(d),cfg,"r")
            buf=StringIO()
            with contextlib.redirect_stdout(buf):
                slurm(run,cfg,True)
            out=buf.getvalue()
            self.assertIn("cg-controller",out)
            self.assertNotIn("cg-preprocess",out)
            self.assertNotIn("cg-kraken_db_setup",out)

    def test_preprocess_cannot_build_missing_kraken_db(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"runs"/"r"; (run/"provenance").mkdir(parents=True)
            cfg={"TAXONOMY_MODE":"kraken2","KRAKEN2_DB":"","KRAKEN2_AUTO_DOWNLOAD":"true","KRAKEN2_DATABASE_SIZE":"standard-8"}
            rows=[{"isolate_id":"i1","group_id":"g","grouping_source":"manifest_group_id","R1":"r1","R2":"r2"}]
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","grouping_source","R1","R2"],rows)
            with self.assertRaisesRegex(SystemExit,"kraken_db_setup"):
                ensure_kraken2_db(run,cfg,rows)

    def test_available_slots_respects_headroom(self):
        self.assertEqual(available_slots(2000,10,0),1990)
        self.assertEqual(available_slots(2000,10,1800),190)
        self.assertEqual(available_slots(2000,10,1990),0)
        self.assertEqual(available_slots(2000,10,2000),0)

    def test_array_task_count(self):
        self.assertEqual(array_task_count("0-499%100"),500)
        self.assertEqual(array_task_count("1,3,7%100"),3)
        self.assertEqual(array_task_count(None),1)

    def test_user_job_count_counts_squeue_rows_for_all_user_jobs(self):
        done=subprocess.CompletedProcess(["squeue"],0,stdout="1\n2\n3\n",stderr="")
        with patch("cleangene.slurm.subprocess.run",return_value=done):
            self.assertEqual(user_job_count("andriy"),3)

    def test_qos_retry_recalculates_and_resubmits(self):
        err=RuntimeError("sbatch failed\nstderr: QOSMaxSubmitJobPerUserLimit")
        with patch("cleangene.slurm.user_job_count",side_effect=[1990,1800,1800]), \
             patch("cleangene.slurm.time.sleep"), \
             patch("cleangene.slurm.submit",side_effect=[err,"123"]):
            with contextlib.redirect_stdout(StringIO()):
                self.assertEqual(submit_with_qos_retry(["sbatch"],{"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0"},1,"x"),"123")

    def test_controller_does_not_oversubmit_when_partially_full(self):
        seen=[]
        cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0","SLURM_MAX_PARALLEL":"100","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
        with tempfile.TemporaryDirectory() as d, \
             patch("cleangene.workers.user_job_count",return_value=1800), \
             patch("cleangene.workers._wait_jobs"), \
             patch("cleangene.workers.submit_with_qos_retry",side_effect=lambda cmd,cfg,task_count,label:(seen.append((cmd,task_count)) or "1")):
            with contextlib.redirect_stdout(StringIO()):
                _run_array_stage(Path(d),cfg,"preprocess",500,"1","1G","1:00:00","x")
        array=seen[0][0][seen[0][0].index("--array")+1]
        self.assertEqual(array,"0-189%100")
        self.assertEqual(seen[0][1],190)

    def test_wait_jobs_sleep_branch_has_time_import(self):
        cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0"}
        with patch("cleangene.workers.job_active",side_effect=[{"1"},set()]), \
             patch("cleangene.workers.user_job_count",return_value=0), \
             patch("cleangene.workers.assert_jobs_succeeded"), \
             patch("cleangene.workers.time.sleep") as sleep, \
             contextlib.redirect_stdout(StringIO()):
            _wait_jobs(["1"],cfg,"x","0/1 complete")
        sleep.assert_called_once_with(0)

    def test_panaroo_creates_nested_regrouped_directory(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; group="Streptococcus gallolyticus"; old="__kraken_pending__"
            pan=Path(d)/"external"; pan.mkdir()
            (pan/"gene_presence_absence.csv").write_text("Gene,iso1,iso2\ng1,l1,l2\n")
            (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            atomic_json(run/"provenance"/"resolved_config.json",{"PANAROO_CPUS":"1","PANAROO_CLEAN_MODE":"strict"})
            rows=[
                {"isolate_id":"iso1","group_id":group,"R1":"r1","R2":"r2","pangenome_dir":str(pan)},
                {"isolate_id":"iso2","group_id":group,"R1":"r1","R2":"r2","pangenome_dir":str(pan)},
            ]
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2","pangenome_dir"],rows)
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[[group,2,"small"]])
            for iso in ("iso1","iso2"):
                qc=run/"results"/"groups"/old/"01_isolates"/iso/"qc.tsv"
                write_tsv(qc,["isolate_id","group_id","excluded","assembly","gff","R1","R2"],[{"isolate_id":iso,"group_id":old,"excluded":0,"assembly":"","gff":"","R1":"r1","R2":"r2"}])
            panaroo(run,0)
            self.assertTrue((run/"results"/"groups"/"Streptococcus_gallolyticus"/"logs").is_dir())

    def test_resume_skips_completed_preprocess_and_resolve(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg={"TAXONOMY_MODE":"off","SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0","SLURM_MAX_PARALLEL":"100","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_ACCOUNT":"","SLURM_PARTITION":"","SLURM_CPUS":"1","SLURM_MEM":"1G","SLURM_TIME":"1:00:00","GROUP_ORCHESTRATOR_CPUS":"1","GROUP_ORCHESTRATOR_MEM":"1G","GROUP_ORCHESTRATOR_TIME":"1:00:00","SUMMARY_CPUS":"1","SUMMARY_MEM":"1G","SUMMARY_TIME":"1:00:00","VALIDATION_CPUS":"1","VALIDATION_MEM":"1G","VALIDATION_TIME":"1:00:00","PANAROO_CPUS":"1","PANAROO_MEM":"1G","PANAROO_TIME":"1:00:00"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id"],[["iso1","g"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["g","iso1"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g",1,"small"]])
            atomic_json(run/"state"/"preprocess"/"iso1.done.json",{})
            atomic_json(run/"state"/"resolve_groups.done.json",{})
            with patch("cleangene.workers._run_index_stage") as arrays, patch("cleangene.workers._run_single_job") as singles, patch("cleangene.workers.controller_downstream"):
                slurm_controller(run)
            self.assertFalse(any("preprocess" in str(c) for c in arrays.call_args_list))
            self.assertFalse(any("resolve_groups" in str(c) for c in singles.call_args_list))

    def test_resume_reruns_failed_panaroo_only(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0","SLURM_MAX_PARALLEL":"100","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_ACCOUNT":"","SLURM_PARTITION":"","PANAROO_CPUS":"1","PANAROO_MEM":"1G","PANAROO_TIME":"1:00:00","PANAROO_SMALL_CPUS":"1","PANAROO_SMALL_MEM":"1G","PANAROO_SMALL_TIME":"1:00:00","SUMMARY_CPUS":"1","SUMMARY_MEM":"1G","SUMMARY_TIME":"1:00:00","VALIDATION_CPUS":"1","VALIDATION_MEM":"1G","VALIDATION_TIME":"1:00:00"}
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id"],[["iso1","done"],["iso2","failed"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["done","iso1"],["failed","iso2"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["done",2,"small"],["failed",2,"small"]])
            atomic_json(run/"state"/"panaroo"/"done.done.json",{})
            seen=[]
            with patch("cleangene.workers._run_index_stage",side_effect=lambda run,cfg,stage,indices,*args:(seen.append((stage,indices)) or [])), patch("cleangene.workers._run_single_job"):
                controller_downstream(run)
            self.assertIn(("panaroo",[1]),seen)
            self.assertNotIn(("panaroo",[0]),seen)

if __name__ == "__main__": unittest.main()
