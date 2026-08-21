import os, tempfile, unittest
from pathlib import Path
from io import StringIO
import contextlib
from cleangene.manifest import groups, load_manifest
from cleangene.evidence import classify_gene_evidence, fixed_coordinate_identity
from cleangene.fasta import assembly_metrics
from cleangene.pangenome import normalize_panaroo, present, select_rows
from cleangene.slurm import array_task_count, available_slots, submit_with_qos_retry, user_job_count, sbatch_cmd, submit
from cleangene.util import atomic_json, read_tsv, write_tsv
from cleangene.workers import _RollingScheduler, _controller_pipeline, _preprocess_scratch, _wait_jobs, controller_downstream, ensure_kraken2_db, kraken_db_for_worker, manifest_pangenome_dir, manifest_row_for_task, panaroo, parse_kraken_report, plot_group, preprocess, reduce_group, slurm_controller, task_row
from cleangene.cli import invalidate_legacy_identity_metrics, make_run, refresh_resume_config, slurm
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

    def test_validation_classification_preserves_unresolved_identity(self):
        d=classify_gene_evidence(mapped_reads=274,breadth=0.998255,mean_depth=59.8726,identity=None,min_breadth=0.90,min_depth=5,min_identity=0.95)
        self.assertEqual(d["validation_state"],"identity_unresolved")
        self.assertEqual(d["validated_call"],"")
        self.assertEqual(d["final_call_source"],"initial_call_unresolved")

    def test_validation_classification_states(self):
        args=dict(min_breadth=0.90,min_depth=5,min_identity=0.95)
        self.assertEqual(classify_gene_evidence(mapped_reads=0,breadth=0,mean_depth=0,identity=None,**args)["validation_state"],"not_detected")
        self.assertEqual(classify_gene_evidence(mapped_reads=10,breadth=0.5,mean_depth=10,identity=1.0,**args)["validation_state"],"partial_coverage")
        self.assertEqual(classify_gene_evidence(mapped_reads=10,breadth=1.0,mean_depth=2,identity=1.0,**args)["validation_state"],"low_depth")
        self.assertEqual(classify_gene_evidence(mapped_reads=10,breadth=1.0,mean_depth=10,identity=0.8,**args)["validation_state"],"divergent")
        self.assertEqual(classify_gene_evidence(mapped_reads=10,breadth=1.0,mean_depth=10,identity=0.99,**args)["validation_state"],"confirmed_present")

    def test_fixed_coordinate_identity_avoids_false_zero_for_matching_names(self):
        self.assertEqual(fixed_coordinate_identity("ACGT","ACGT")["identity"],1.0)
        self.assertEqual(fixed_coordinate_identity("ACGT","ACGA")["identity"],0.75)
        self.assertIsNone(fixed_coordinate_identity("NNNN","NNNN"))

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
        cmd=sbatch_cmd(name="x",wrap="echo ok",cpus="1",mem="1G",time="00:05:00",array="0-499%100")
        self.assertIn("--array",cmd)
        self.assertEqual(cmd[cmd.index("--array")+1],"0-499%100")

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

    def test_preprocess_threads_kraken_and_suppresses_per_read_output(self):
        with tempfile.TemporaryDirectory() as d:
            run_dir=Path(d)/"run"; db=Path(d)/"db"; db.mkdir(); (db/"hash.k2d").write_text("db")
            r1=Path(d)/"r1.fq"; r2=Path(d)/"r2.fq"; r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n")
            cfg={"TAXONOMY_MODE":"kraken2","KRAKEN2_DB":str(db),"KRAKEN2_DB_ACCESS":"direct","KRAKEN2_KEEP_CLASSIFICATIONS":"false","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"true","READ_TRIMMING_MODE":"off","CPUS":"8"}
            rows=[["i1","g","manifest_organism","Species",str(r1),str(r2),"/prebuilt"]]
            atomic_json(run_dir/"provenance"/"resolved_config.json",cfg); write_tsv(run_dir/"provenance"/"manifest.tsv",["isolate_id","group_id","grouping_source","organism","R1","R2","pangenome_dir"],rows); write_tsv(run_dir/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["g","i1"]])
            commands=[]
            def fake_run(command,**kwargs):
                commands.append(command)
                if command[0]=="kraken2": Path(command[command.index("--report")+1]).write_text("100.00\t1\t1\tS\t1\tSpecies\n")
            with patch("cleangene.workers.run",side_effect=fake_run), patch.dict(os.environ,{"SLURM_TMPDIR":str(Path(d)/"scratch")},clear=False): preprocess(run_dir,0)
            kraken=next(c for c in commands if c[0]=="kraken2")
            self.assertEqual(kraken[kraken.index("--threads")+1],"8")
            self.assertEqual(kraken[kraken.index("--output")+1],"/dev/null")
            shared=run_dir/"results"/"groups"/"g"/"01_isolates"/"i1"
            self.assertFalse((shared/"kraken2.output.tsv").exists()); self.assertTrue((shared/"kraken2.report.tsv").is_file()); self.assertTrue((shared/"qc.tsv").is_file())

    def test_kraken_database_stages_once_in_node_cache(self):
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/"source"; cache=Path(d)/"cache"; source.mkdir(); (source/"hash.k2d").write_text("hash"); (source/"opts.k2d").write_text("opts")
            cfg={"KRAKEN2_DB_ACCESS":"copy","KRAKEN2_NODE_CACHE_DIR":str(cache),"KRAKEN2_NODE_CACHE_MIN_FREE_GB":"0"}
            first,mmap1=kraken_db_for_worker(str(source),cfg); second,mmap2=kraken_db_for_worker(str(source),cfg)
            self.assertEqual(first,second); self.assertTrue(mmap1 and mmap2); self.assertTrue((Path(first)/".cleangene-ready").is_file())

    def test_preprocess_scratch_uses_slurm_tmpdir(self):
        with tempfile.TemporaryDirectory() as d, patch.dict(os.environ,{"SLURM_TMPDIR":d},clear=False):
            scratch=_preprocess_scratch({"PREPROCESS_USE_NODE_LOCAL_SCRATCH":"true","PREPROCESS_SCRATCH_DIR":""},Path(d)/"run","iso")
            self.assertEqual(scratch.parent,Path(d)); self.assertTrue(scratch.is_dir())

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
        with patch("cleangene.slurm.subprocess.run",return_value=done) as command:
            self.assertEqual(user_job_count("andriy"),3)
        self.assertIn("%F|%T",command.call_args.args[0])

    def test_qos_retry_recalculates_and_resubmits(self):
        err=RuntimeError("sbatch failed\nstderr: QOSMaxSubmitJobPerUserLimit")
        with patch("cleangene.slurm.user_job_count",side_effect=[1990,1800,1800]), \
             patch("cleangene.slurm.time.sleep"), \
             patch("cleangene.slurm.submit",side_effect=[err,"123"]):
            with contextlib.redirect_stdout(StringIO()):
                self.assertEqual(submit_with_qos_retry(["sbatch"],{"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0"},1,"x"),"123")

    def test_controller_does_not_oversubmit_when_partially_full(self):
        seen=[]
        cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0","SLURM_MAX_PARALLEL":"400","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_MAX_OUTSTANDING_CHUNKS":"8","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
        with tempfile.TemporaryDirectory() as d, \
             patch("cleangene.workers.submit_with_qos_retry",side_effect=lambda cmd,cfg,task_count,label:(seen.append((cmd,task_count)) or "1")):
            run=Path(d); write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],(("g",f"i{i}") for i in range(500)))
            scheduler=_RollingScheduler(run,cfg); scheduler.snapshot={"total":1800,"jobs":{}}
            with contextlib.redirect_stdout(StringIO()):
                scheduler.submit_ready("preprocess",list(range(500)),"1","1G","1:00:00","x",400)
        array=seen[0][0][seen[0][0].index("--array")+1]
        self.assertEqual(array,"0-189%400")
        self.assertEqual(seen[0][1],190)

    def test_rolling_scheduler_refills_before_prior_array_finishes(self):
        cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_MAX_PARALLEL":"400","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_MAX_OUTSTANDING_CHUNKS":"8","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],(("g",f"i{i}") for i in range(600)))
            scheduler=_RollingScheduler(run,cfg)
            with patch("cleangene.workers.submit_with_qos_retry",side_effect=["10","11"]):
                scheduler.snapshot={"total":0,"jobs":{}}
                self.assertEqual(scheduler.submit_ready("preprocess",list(range(600)),"1","1G","1:00:00","x",400),400)
                for i in range(200): atomic_json(run/"state"/"preprocess"/f"i{i}.done.json",{})
                scheduler.snapshot={"total":200,"jobs":{"10":{"RUNNING":200}}}
                self.assertEqual(scheduler.submit_ready("preprocess",list(range(600)),"1","1G","1:00:00","x",400),200)
            self.assertIn("10",scheduler.active)
            self.assertIn("11",scheduler.active)

    def test_progress_counts_done_markers_not_submitted_boundary(self):
        cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10"}
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],(("g",f"i{i}") for i in range(3))); atomic_json(run/"state"/"preprocess"/"i0.done.json",{})
            scheduler=_RollingScheduler(run,cfg); scheduler.snapshot={"total":2,"jobs":{}}; scheduler.submitted={"preprocess":{0,1,2}}
            output=StringIO()
            with contextlib.redirect_stdout(output): scheduler.progress("preprocess",[0,1,2],"CleanGene preprocess")
            self.assertIn("complete=1/3 submitted=3",output.getvalue())

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
            with patch("cleangene.workers._controller_pipeline") as pipeline, patch("cleangene.workers._run_single_job") as singles:
                slurm_controller(run)
            pipeline.assert_called_once_with(run,True)
            self.assertFalse(any("resolve_groups" in str(c) for c in singles.call_args_list))

    def test_resume_reruns_failed_panaroo_only(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0","SLURM_MAX_PARALLEL":"100","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_ACCOUNT":"","SLURM_PARTITION":"","PANAROO_CPUS":"1","PANAROO_MEM":"1G","PANAROO_TIME":"1:00:00","PANAROO_SMALL_CPUS":"1","PANAROO_SMALL_MEM":"1G","PANAROO_SMALL_TIME":"1:00:00","SUMMARY_CPUS":"1","SUMMARY_MEM":"1G","SUMMARY_TIME":"1:00:00","VALIDATION_CPUS":"1","VALIDATION_MEM":"1G","VALIDATION_TIME":"1:00:00","PLOT_CPUS":"1","PLOT_MEM":"1G","PLOT_TIME":"1:00:00"}
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

    def test_known_groups_overlap_downstream_with_preprocess(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; atomic_json(run/"provenance"/"resolved_config.json",{})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id"],[["small_iso","small"],["large_iso","large"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["small","small_iso"],["large","large_iso"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["small",1,"small"],["large",2001,"large"]])
            atomic_json(run/"state"/"preprocess"/"small_iso.done.json",{})
            class ImmediateScheduler:
                def __init__(self,*args): self.active={}; self.calls=[]
                def refresh(self): pass
                def active_indices(self,stage): return set()
                def submit_ready(self,stage,indices,*args):
                    self.calls.append((stage,list(indices)))
                    kind="isolate" if stage in {"preprocess","validate"} else "group"; rows=read_tsv(run/"state"/f"{kind}_tasks.tsv")
                    for i in indices:
                        name=rows[i]["isolate_id" if kind=="isolate" else "group_id"]; atomic_json(run/"state"/stage/f"{name}.done.json",{})
                    return len(indices)
                def progress(self,*args): pass
                def wait_tick(self): pass
            scheduler=ImmediateScheduler()
            with patch("cleangene.workers._RollingScheduler",return_value=scheduler), patch("cleangene.workers._run_single_job",return_value="summary"):
                _controller_pipeline(run,True)
            panaroo_calls=[indices for stage,indices in scheduler.calls if stage=="panaroo"]
            self.assertEqual(panaroo_calls[0],[0])
            self.assertEqual(panaroo_calls[1],[1])
            self.assertLess(scheduler.calls.index(("panaroo",[0])),scheduler.calls.index(("preprocess",[0,1])))

    def test_resume_config_updates_execution_values_and_preserves_db(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; config=Path(d)/"arc.env"; atomic_json(run/"provenance"/"resolved_config.json",{"KRAKEN2_DB":"/resolved/db","SLURM_MAX_PARALLEL":"100"})
            config.write_text('SLURM_MAX_PARALLEL="400"\nSLURM_PREPROCESS_MAX_INFLIGHT="400"\nKRAKEN2_DB=""\n')
            updated=refresh_resume_config(run,config)
            self.assertEqual(updated["SLURM_MAX_PARALLEL"],"400"); self.assertEqual(updated["KRAKEN2_DB"],"/resolved/db")
            self.assertTrue((run/"provenance"/"resolved_config.pre_resume.json").is_file())

    def test_reduce_preserves_initial_call_for_identity_unresolved(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; group="g"; iso="iso1"; other="iso2"; root=run/"results"/"groups"/group
            (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            atomic_json(run/"provenance"/"resolved_config.json",{})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[[iso,group,"r1","r2"],[other,group,"r1","r2"]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[[group,2,"small"]])
            for sample in (iso,other):
                write_tsv(root/"01_isolates"/sample/"qc.tsv",["isolate_id","group_id","excluded","assembly","gff","R1","R2"],[[sample,group,0,"","", "r1","r2"]])
            write_tsv(root/"02_pangenome"/"initial_calls"/"gene_presence_absence.binary.tsv",["Gene",iso,other],[["rpoE",1,1]])
            write_tsv(root/"03_read_validation"/"evidence"/iso/"metrics.tsv",["reference_id","Gene","validated_call","validation_state","decision_reason","final_call_source","breadth","mean_depth","identity","mapped_reads"],[["rpoE","rpoE","","identity_unresolved","identity could not be measured","initial_call_unresolved",0.99,50,"NA",274]])
            reduce_group(run,0)
            rows=read_tsv(root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv")
            self.assertEqual(rows[0][iso],"1")

    def test_plot_group_writes_done_marker(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; group="g"; root=run/"results"/"groups"/group
            (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            atomic_json(run/"provenance"/"resolved_config.json",{"PLOT_MAX_CLUSTER_ISOLATES":"2000"})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id"],[["i",group]])
            write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[[group,1,"small"]])
            write_tsv(root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv",["Gene","i"],[["g1",1]])
            with patch("cleangene.workers.plot_presence_absence") as plot:
                plot_group(run,0)
            plot.assert_called_once()
            self.assertTrue((run/"state"/"plot"/"g.done.json").is_file())

    def test_resume_invalidates_legacy_zero_identity_metrics(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; g="g"; iso="iso"; metrics=run/"results"/"groups"/g/"03_read_validation"/"evidence"/iso/"metrics.tsv"
            (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            write_tsv(metrics,["reference_id","Gene","validated_call","validation_state","breadth","mean_depth","identity","aligned_positions","mapped_reads"],[["rpoE","rpoE",0,"partial_or_divergent",0.99,50,0,0,274]])
            for p in (run/"state"/"validate"/"iso.done.json",run/"state"/"reduce"/"g.done.json",run/"state"/"plot"/"g.done.json",run/"state"/"summary.done.json"):
                atomic_json(p,{})
            n=invalidate_legacy_identity_metrics(run,{"READ_VALIDATION_MIN_BREADTH":"0.90","READ_VALIDATION_MIN_MEAN_DEPTH":"5"})
            self.assertEqual(n,1)
            self.assertTrue(metrics.with_name("metrics.pre_identity_fix.tsv").is_file())
            self.assertFalse((run/"state"/"validate"/"iso.done.json").exists())
            self.assertFalse((run/"state"/"reduce"/"g.done.json").exists())

if __name__ == "__main__": unittest.main()
