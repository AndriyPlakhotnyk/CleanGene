import gzip, os, tempfile, unittest
from pathlib import Path
from io import StringIO
import contextlib
from cleangene.manifest import groups, load_manifest
from cleangene.evidence import classify_gene_evidence, fixed_coordinate_identity
from cleangene.fasta import assembly_metrics
from cleangene.pangenome import normalize_panaroo, present, select_rows
from cleangene.slurm import array_task_count, available_slots, submit_with_qos_retry, user_job_count, user_queue_snapshot, sbatch_cmd, submit
from cleangene.task_store import build_isolate_task_store, load_isolate_task, migrate_isolate_task_store
from cleangene.util import atomic_json, read_tsv, write_tsv
from cleangene.workers import _RollingScheduler, _controller_cmd, _controller_pipeline, _index_done, _preprocess_scratch, _wait_jobs, build_organism_results_index, cleanup_trimmed_fastqs, compress_completed_outputs, controller_downstream, ensure_kraken2_db, kraken_db_for_worker, manifest_pangenome_dir, manifest_row_for_task, panaroo, parse_kraken_report, plot_group, prepare_read_inputs, preprocess, reduce_group, slurm_controller, task_row
from cleangene.cli import apply_cli_overrides, exclude_command, invalidate_legacy_identity_metrics, local, make_run, refresh_resume_config, slurm
from unittest.mock import patch
import subprocess

class CleanGeneCoreTests(unittest.TestCase):
    def test_organism_results_index_links_complete_isolate_directories(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; source=run/"results"/"groups"/"pending"/"01_isolates"/"BI_06_0500"
            (source/"assembly").mkdir(parents=True); (source/"annotation").mkdir(); (source/"logs").mkdir()
            (source/"assembly"/"contigs.fa").write_text(">c\nACGT\n"); (source/"annotation"/"BI_06_0500.gff").write_text("##gff-version 3\n"); (source/"logs"/"prokka.stderr").write_text("")
            atomic_json(run/"provenance"/"resolved_config.json",{})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","organism"],[["BI_06_0500","pending","Expected species"]])
            write_tsv(source/"qc.tsv",["isolate_id","top_species"],[["BI_06_0500","Identified species"]])
            result=build_organism_results_index(run); link=run/"results"/"organisms"/"Identified_species"/"BI_06_0500"
            self.assertEqual(result,{"organisms":1,"isolates":1}); self.assertTrue(link.is_symlink()); self.assertEqual(link.resolve(),source.resolve())
            self.assertEqual((link/"assembly"/"contigs.fa").read_text(),">c\nACGT\n"); self.assertTrue((link/"annotation"/"BI_06_0500.gff").is_file()); self.assertTrue((link/"logs"/"prokka.stderr").is_file())
            index=read_tsv(run/"results"/"cohort"/"organism_isolate_index.tsv"); self.assertEqual(index[0]["organism"],"Identified species"); self.assertEqual(index[0]["isolate_id"],"BI_06_0500")

    def test_organism_results_index_refreshes_stale_links_and_uses_manifest_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; source=run/"results"/"groups"/"g"/"01_isolates"/"BI_1"; source.mkdir(parents=True)
            atomic_json(run/"provenance"/"resolved_config.json",{}); write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","organism"],[["BI_1","g","Manifest species"]]); write_tsv(source/"qc.tsv",["isolate_id","top_species"],[["BI_1",""]])
            stale=run/"results"/"organisms"/"Old_species"/"BI_1"; stale.parent.mkdir(parents=True); stale.symlink_to(source)
            build_organism_results_index(run)
            self.assertFalse(stale.is_symlink()); self.assertTrue((run/"results"/"organisms"/"Manifest_species"/"BI_1").is_symlink()); self.assertFalse(stale.parent.exists())

    def test_cleanup_replaces_trimmed_fastqs_with_original_links(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; group="g"; iso="iso1"; original=Path(d)/"input"; original.mkdir()
            r1=original/"iso1_R1.fastq.gz"; r2=original/"iso1_R2.fastq.gz"; r1.write_bytes(b"original-r1"); r2.write_bytes(b"original-r2")
            atomic_json(run/"provenance"/"resolved_config.json",{})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[[iso,group,r1,r2]])
            reads=run/"results"/"groups"/group/"01_isolates"/iso/"reads"; reads.mkdir(parents=True)
            tr1=reads/"trimmed_R1.fastq.gz"; tr2=reads/"trimmed_R2.fastq.gz"; tr1.write_bytes(b"trimmed-r1"); tr2.write_bytes(b"trimmed-r2")
            write_tsv(reads.parent/"qc.tsv",["isolate_id","group_id","excluded","R1","R2"],[[iso,group,0,tr1,tr2]])
            preview=cleanup_trimmed_fastqs(run,dry_run=True)
            self.assertEqual(preview["counts"],{"would_link":1}); self.assertFalse(tr1.is_symlink())
            result=cleanup_trimmed_fastqs(run)
            self.assertEqual(result["counts"],{"linked":1}); self.assertTrue(tr1.is_symlink()); self.assertTrue(tr2.is_symlink())
            self.assertEqual(tr1.read_bytes(),r1.read_bytes()); self.assertEqual(tr2.read_bytes(),r2.read_bytes())
            self.assertTrue((run/"results"/"cohort"/"fastq_cleanup.tsv").is_file())

    def test_final_storage_sweep_compresses_previously_completed_preprocess(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; iso_dir=run/"results"/"groups"/"g"/"01_isolates"/"iso1"
            assembly=iso_dir/"assembly"; annotation=iso_dir/"annotation"; assembly.mkdir(parents=True); annotation.mkdir()
            contigs=assembly/"contigs.fa"; contigs.write_text(">c\nACGT\n"); (assembly/"spades.fasta").write_text(">s\nACGT\n"); (assembly/"spades.gfa").write_text("H\tVN:Z:1.0\n")
            gff=annotation/"iso1.gff"; gff.write_text("##gff-version 3\n"); (annotation/"iso1.gbk").write_text("LOCUS\n"); (annotation/"iso1.sqn").write_text("sqn\n")
            atomic_json(run/"provenance"/"resolved_config.json",{"COMPRESS_ASSEMBLY_OUTPUTS":"intermediates","COMPRESS_ANNOTATION_OUTPUTS":"nonessential"})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["iso1","g","r1","r2"]])
            write_tsv(iso_dir/"qc.tsv",["isolate_id","group_id","assembly","gff"],[["iso1","g",contigs,gff]])
            result=compress_completed_outputs(run)
            self.assertEqual(result["files_compressed"],4)
            self.assertTrue(contigs.is_file())
            self.assertTrue((assembly/"spades.fasta.gz").is_file()); self.assertTrue((assembly/"spades.gfa.gz").is_file())
            self.assertTrue(gff.is_file()); self.assertTrue((annotation/"iso1.gbk.gz").is_file()); self.assertTrue((annotation/"iso1.sqn.gz").is_file())
            self.assertEqual(read_tsv(iso_dir/"qc.tsv")[0]["assembly"],str(contigs))
            self.assertTrue((run/"results"/"cohort"/"storage_cleanup.tsv").is_file())

    def test_cleanup_flag_sets_final_summary_cleanup(self):
        args=type("Args",(),{"cleanup_trimmed_fastq":True})()
        self.assertEqual(apply_cli_overrides({},args)["CLEANUP_TRIMMED_FASTQ"],"true")

    def test_exclude_preserves_manifest_row_and_skips_unfinished_preprocess(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"; r1.write_text("reads1"); r2.write_text("reads2")
            manifest=root/"manifest.tsv"; manifest.write_text(f"isolate_id\tR1\tR2\niso1\t{r1}\t{r2}\n")
            run=make_run(manifest,root,{"TAXONOMY_MODE":"off","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false"},"r")
            args=type("Args",(),{"run_dir":run,"samples":["iso1"],"samples_file":None})()
            exclude_command(args)
            rows=read_tsv(run/"provenance"/"manifest.tsv"); self.assertEqual(len(rows),1); self.assertEqual(rows[0]["user_excluded"],"true"); self.assertEqual(rows[0]["grouping_source"],"kraken_pending")
            with patch("cleangene.workers.run",side_effect=AssertionError("excluded isolate invoked an external tool")): preprocess(run,0)
            qc=read_tsv(run/"results"/"groups"/"__kraken_pending__"/"01_isolates"/"iso1"/"qc.tsv")[0]
            self.assertEqual(qc["excluded"],"1"); self.assertEqual(qc["reason"],"user_excluded")
            from cleangene.workers import resolve_groups
            resolve_groups(run,None); resolved=read_tsv(run/"provenance"/"manifest.tsv")[0]; self.assertEqual(resolved["group_id"],"__user_excluded__")

    def test_direct_spades_uses_original_symlinked_reads(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"; sequence="A"*120; quality="I"*120; r1.write_text(f"@r\n{sequence}\n+\n{quality}\n"); r2.write_text(f"@r\n{sequence}\n+\n{quality}\n")
            manifest=root/"manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            cfg={"TAXONOMY_MODE":"off","ASSEMBLER":"spades","CHECKM2_MODE":"off","QC_MIN_N50_PASS":"0","QC_MIN_N50_FAIL":"0","QC_MIN_COVERAGE_PASS":"0","QC_MIN_COVERAGE_FAIL":"0","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false"}; run_dir=make_run(manifest,root,cfg,"r"); commands=[]
            def fake_run(command,**kwargs):
                commands.append(command)
                if command[0]=="spades.py":
                    out=Path(command[command.index("-o")+1]); out.mkdir(parents=True,exist_ok=True); (out/"contigs.fasta").write_text(">c\nACGT\n"); (out/"assembly_graph.gfa").write_text("H\tVN:Z:1.0\n")
                elif command[0]=="prokka":
                    out=Path(command[command.index("--outdir")+1]); prefix=command[command.index("--prefix")+1]; out.mkdir(parents=True,exist_ok=True); (out/f"{prefix}.gff").write_text("##gff-version 3\n")
            with patch("cleangene.workers.run",side_effect=fake_run): preprocess(run_dir,0)
            self.assertEqual(sum(command[0]=="spades.py" for command in commands),1); self.assertFalse(any(command[0]=="shovill" for command in commands))
            spades=next(command for command in commands if command[0]=="spades.py"); self.assertIn("--only-assembler",spades); self.assertTrue(Path(spades[spades.index("-1")+1]).is_symlink()); self.assertTrue(Path(spades[spades.index("-2")+1]).is_symlink())
            qc=read_tsv(run_dir/"results"/"groups"/"g"/"01_isolates"/"iso1"/"qc.tsv")[0]; self.assertEqual(Path(qc["assembly"]).name,"contigs.fasta"); self.assertIn("+symlinked",qc["read_preprocessing"])
            for field in ("PASS/FAIL","Notes","trimmed_read_length","mean_base_quality","sequencing_coverage","checkm2_completeness","checkm2_contamination","qc_profile_source"): self.assertIn(field,qc)

    def test_exclude_filters_an_already_preprocessed_isolate(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); manifest=root/"manifest.tsv"; manifest.write_text("isolate_id\tgroup_id\tR1\tR2\niso1\tg\tr1\tr2\n")
            run=make_run(manifest,root,{},"r"); iso=run/"results"/"groups"/"g"/"01_isolates"/"iso1"; gff=iso/"annotation"/"iso1.gff"; gff.parent.mkdir(parents=True); gff.write_text("##gff-version 3\n")
            write_tsv(iso/"qc.tsv",["isolate_id","group_id","excluded","R1","R2","assembly","gff"],[["iso1","g",0,"r1","r2","",gff]]); atomic_json(run/"state"/"preprocess"/"iso1.done.json",{"excluded":False,"gff":str(gff)})
            exclude_command(type("Args",(),{"run_dir":run,"samples":["iso1"],"samples_file":None})())
            from cleangene.workers import retained_rows
            self.assertEqual(retained_rows(run,"g"),[]); self.assertTrue(gff.is_file())

    def test_preprocess_slurm_logs_use_stage_subdirectory(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; cfg={"SLURM_ACCOUNT":"","SLURM_PARTITION":""}
            command=_controller_cmd(run,cfg,"preprocess","1-2%2","1","1G","01:00:00")
            log=Path(command[command.index("--output")+1]); self.assertEqual(log.parent,run/"logs"/"slurm"/"preprocess"); self.assertTrue(log.parent.is_dir())

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
            gz=Path(d)/"x.fa.gz"
            with gzip.open(gz,"wt") as handle: handle.write(p.read_text())
            self.assertEqual(assembly_metrics(gz),m)

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

    def test_indexed_task_store_loads_isolate_by_offset(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); (run/"provenance").mkdir(parents=True)
            rows=[{"isolate_id":f"iso{i}","group_id":"g","organism":"Species A","R1":f"r{i}_1.fq","R2":f"r{i}_2.fq"} for i in range(3)]
            write_tsv(run/"provenance"/"qc_thresholds.tsv",["isolate_id","qc_profile_source","qc_max_contigs_pass","qc_max_contigs_fail","qc_min_n50_pass","qc_min_n50_fail","qc_min_coverage_pass","qc_min_coverage_fail","qc_min_read_length_pass","qc_min_read_length_fail","qc_min_mean_base_quality_pass","qc_min_mean_base_quality_fail","qc_min_completeness_pass","qc_min_completeness_fail","qc_max_checkm2_contamination_pass","qc_max_checkm2_contamination_fail","qc_max_kraken_contamination_fail"],[[r["isolate_id"],"global",300,1000,25000,5000,20,10,120,"",30,"",90,80,5,10,5] for r in rows])
            self.assertEqual(build_isolate_task_store(run,rows),3)
            record=load_isolate_task(run,2)
            self.assertEqual(record["isolate_id"],"iso2")
            self.assertEqual(record["qc_thresholds_resolved"]["qc_max_contigs_pass"],300.0)

    def test_legacy_run_task_store_migrates_from_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); (run/"provenance").mkdir(parents=True)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["iso1","g","r1","r2"]])
            self.assertEqual(migrate_isolate_task_store(run),1)
            self.assertEqual(load_isolate_task(run,0)["isolate_id"],"iso1")

    def test_reads_processed_auto_skips_fastp_but_always_runs(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fq"; r2=root/"r2.fq"; logs=root/"logs"; out=root/"out"
            r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n"); logs.mkdir()
            row={"isolate_id":"iso1","R1":str(r1),"R2":str(r2),"reads_processed":"true"}
            with patch("cleangene.workers.command_exists",return_value=True), patch("cleangene.workers.run") as tool_run:
                _,_,_,trimmed,decision=prepare_read_inputs(row,out,logs,{"READ_TRIMMING_MODE":"auto","SKIP_TRIM":"false","ASSEMBLER":"shovill"})
                self.assertEqual(trimmed,0); self.assertEqual(decision,"skipped_reads_processed_true"); tool_run.assert_not_called()
            with patch("cleangene.workers.command_exists",return_value=True), patch("cleangene.workers.run") as tool_run:
                _,_,_,trimmed,decision=prepare_read_inputs(row,root/"out2",logs,{"READ_TRIMMING_MODE":"always","SKIP_TRIM":"false","ASSEMBLER":"shovill","CPUS":"1"})
                self.assertEqual(trimmed,1); self.assertEqual(decision,"ran_fastp_always"); tool_run.assert_called_once()

    def test_qc_only_skip_shovill_uses_fastq_symlinks(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); manifest=root/"manifest.tsv"
            rows=[]
            for iso in ("iso1","iso2"):
                r1=root/f"{iso}_R1.fastq.gz"; r2=root/f"{iso}_R2.fastq.gz"
                r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n")
                rows.append(f"{iso}\tg\t{r1}\t{r2}")
            manifest.write_text("isolate_id\tgroup_id\tR1\tR2\n" + "\n".join(rows) + "\n")
            cfg={"TAXONOMY_MODE":"off","READ_TRIMMING_MODE":"auto","SKIP_TRIM":"true","SKIP_SHOVILL":"true","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
            run=make_run(manifest,root,cfg,"r")
            local(run)
            qc=read_tsv(run/"results"/"cohort"/"isolate_qc.tsv")
            self.assertEqual(len(qc),2)
            for row in qc:
                self.assertEqual(row["reason"],"shovill_skipped")
                self.assertEqual(row["adapter_trimmed"],"0")
                self.assertIn("+symlinked",row["read_preprocessing"])
                self.assertTrue(Path(row["R1"]).is_symlink())
                self.assertTrue(Path(row["R2"]).is_symlink())
                self.assertEqual(row["assembly"],"")
                self.assertEqual(row["gff"],"")
                for field in ("PASS/FAIL","Notes","trimmed_read_length","mean_base_quality","sequencing_coverage","checkm2_completeness","checkm2_contamination","qc_profile_source"): self.assertIn(field,row)
            self.assertTrue((run/"state"/"summary.done.json").is_file())

    def test_assembly_qc_failure_skips_prokka(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); sequence="A"*120; quality="I"*120; r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"
            r1.write_text(f"@r\n{sequence}\n+\n{quality}\n"); r2.write_text(f"@r\n{sequence}\n+\n{quality}\n")
            manifest=root/"manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            run_dir=make_run(manifest,root,{"TAXONOMY_MODE":"off","READ_TRIMMING_MODE":"off","CHECKM2_MODE":"off","QC_MIN_COVERAGE_PASS":"0","QC_MIN_COVERAGE_FAIL":"0","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false"},"r"); commands=[]
            def fake_run(command,**kwargs):
                commands.append(command)
                if command[0]=="shovill":
                    out=Path(command[command.index("--outdir")+1]); out.mkdir(parents=True,exist_ok=True); (out/"contigs.fa").write_text(">c\nACGT\n")
            with patch("cleangene.workers.run",side_effect=fake_run): preprocess(run_dir,0)
            self.assertFalse(any(command[0]=="prokka" for command in commands))
            qc=read_tsv(run_dir/"results"/"groups"/"g"/"01_isolates"/"iso1"/"qc.tsv")[0]
            self.assertEqual(qc["PASS/FAIL"],"FAIL"); self.assertIn("n50_low",qc["reason"]); self.assertIn("gff_missing",qc["reason"])

    def test_prokka_failure_marks_isolate_fail(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); sequence="A"*120; quality="I"*120; r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"
            r1.write_text(f"@r\n{sequence}\n+\n{quality}\n"); r2.write_text(f"@r\n{sequence}\n+\n{quality}\n")
            manifest=root/"manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            cfg={"TAXONOMY_MODE":"off","READ_TRIMMING_MODE":"off","CHECKM2_MODE":"off","QC_MIN_N50_PASS":"0","QC_MIN_N50_FAIL":"0","QC_MIN_COVERAGE_PASS":"0","QC_MIN_COVERAGE_FAIL":"0","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false"}; run_dir=make_run(manifest,root,cfg,"r")
            def fake_run(command,**kwargs):
                if command[0]=="shovill":
                    out=Path(command[command.index("--outdir")+1]); out.mkdir(parents=True,exist_ok=True); (out/"contigs.fa").write_text(">c\nACGT\n")
                elif command[0]=="prokka": raise subprocess.CalledProcessError(2,command)
            with patch("cleangene.workers.run",side_effect=fake_run): preprocess(run_dir,0)
            qc=read_tsv(run_dir/"results"/"groups"/"g"/"01_isolates"/"iso1"/"qc.tsv")[0]
            self.assertEqual(qc["PASS/FAIL"],"FAIL"); self.assertEqual(qc["excluded"],"1"); self.assertIn("prokka_failed",qc["reason"])

    def test_compress_assembly_outputs_updates_qc_path(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"
            r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n")
            manifest=root/"manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            cfg={"TAXONOMY_MODE":"off","READ_TRIMMING_MODE":"off","COMPRESS_ASSEMBLY_OUTPUTS":"all","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
            run=make_run(manifest,root,cfg,"r")
            def fake_run(command,**kwargs):
                if command[0]=="shovill":
                    out=Path(command[command.index("--outdir")+1]); out.mkdir(parents=True,exist_ok=True)
                    (out/"contigs.fa").write_text(">c1\nACGT\n")
                    (out/"spades.fasta").write_text(">s1\nACGT\n")
                    (out/"spades.gfa").write_text("H\tVN:Z:1.0\n")
                elif command[0]=="prokka":
                    out=Path(command[command.index("--outdir")+1]); prefix=command[command.index("--prefix")+1]; out.mkdir(parents=True,exist_ok=True)
                    (out/f"{prefix}.gff").write_text("##gff-version 3\n")
            with patch("cleangene.workers.run",side_effect=fake_run): preprocess(run,0)
            qc=read_tsv(run/"results"/"groups"/"g"/"01_isolates"/"iso1"/"qc.tsv")[0]
            assembly=Path(qc["assembly"])
            self.assertEqual(assembly.name,"contigs.fa.gz")
            self.assertTrue(assembly.is_file())
            self.assertFalse((assembly.parent/"contigs.fa").exists())
            self.assertTrue((assembly.parent/"spades.fasta.gz").is_file())
            self.assertTrue((assembly.parent/"spades.gfa.gz").is_file())
            self.assertEqual(assembly_metrics(assembly)["assembly_length"],4)

    def test_compress_annotation_outputs_keeps_gff_plain(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"
            r1.write_text("@r\nA\n+\n!\n"); r2.write_text("@r\nT\n+\n!\n")
            manifest=root/"manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            cfg={"TAXONOMY_MODE":"off","READ_TRIMMING_MODE":"off","CHECKM2_MODE":"off","QC_MIN_N50_PASS":"0","QC_MIN_N50_FAIL":"0","QC_MIN_COVERAGE_PASS":"0","QC_MIN_COVERAGE_FAIL":"0","COMPRESS_ANNOTATION_OUTPUTS":"nonessential","PREPROCESS_USE_NODE_LOCAL_SCRATCH":"false","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
            run=make_run(manifest,root,cfg,"r")
            def fake_run(command,**kwargs):
                if command[0]=="shovill":
                    out=Path(command[command.index("--outdir")+1]); out.mkdir(parents=True,exist_ok=True)
                    (out/"contigs.fa").write_text(">c1\nACGT\n")
                elif command[0]=="prokka":
                    out=Path(command[command.index("--outdir")+1]); prefix=command[command.index("--prefix")+1]; out.mkdir(parents=True,exist_ok=True)
                    for ext in ("gff","sqn","gbk","err","ffn","fna"):
                        (out/f"{prefix}.{ext}").write_text(f"{ext}\n")
            with patch("cleangene.workers.run",side_effect=fake_run): preprocess(run,0)
            ann=run/"results"/"groups"/"g"/"01_isolates"/"iso1"/"annotation"
            qc=read_tsv(ann.parent/"qc.tsv")[0]
            self.assertEqual(Path(qc["gff"]).name,"iso1.gff")
            self.assertTrue((ann/"iso1.gff").is_file())
            for ext in ("sqn","gbk","err","ffn","fna"):
                self.assertFalse((ann/f"iso1.{ext}").exists())
                self.assertTrue((ann/f"iso1.{ext}.gz").is_file())

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
        self.assertIn("%F|%K|%T|%j|%C|%o",command.call_args.args[0])

    def test_squeue_snapshot_parses_array_elements_and_states(self):
        stdout="123|7|R|cg-preprocess|8|python -m cleangene _worker\n123|8|PENDING|cg-preprocess|8|python -m cleangene _worker\n999|N/A|RUNNING|unrelated|1|sleep 10\n"
        done=subprocess.CompletedProcess(["squeue"],0,stdout=stdout,stderr="")
        with patch("cleangene.slurm.subprocess.run",return_value=done): snapshot=user_queue_snapshot("andriy")
        self.assertEqual(snapshot["total"],3)
        self.assertEqual(snapshot["jobs"]["123"]["RUNNING"],1)
        self.assertEqual(snapshot["jobs"]["123"]["PENDING"],1)
        self.assertEqual(snapshot["entries"][0]["task_id"],"7")
        self.assertEqual(snapshot["entries"][0]["cpus"],8)

    def test_resumed_controller_adopts_live_run_arrays(self):
        cfg={"SLURM_USER_JOB_LIMIT":"2000","SLURM_JOB_HEADROOM":"10","SLURM_POLL_SECONDS":"0","SLURM_MAX_PARALLEL":"100","SLURM_ARRAY_CHUNK_SIZE":"500","SLURM_MAX_OUTSTANDING_CHUNKS":"8","SLURM_ACCOUNT":"","SLURM_PARTITION":""}
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],(("g",f"i{i}") for i in range(3)))
            snapshot={"total":3,"jobs":{"123":{"RUNNING":1,"PENDING":1},"999":{"RUNNING":1}},"entries":[
                {"job_id":"123","task_id":"1","state":"RUNNING","name":"cg-preprocess","command":f"python -m cleangene _worker --stage preprocess --run-dir {run} --index 1"},
                {"job_id":"123","task_id":"2","state":"PENDING","name":"cg-preprocess","command":f"python -m cleangene _worker --stage preprocess --run-dir {run} --index 2"},
                {"job_id":"999","task_id":"0","state":"RUNNING","name":"cg-preprocess","command":"python -m cleangene _worker --stage preprocess --run-dir /another/run --index 0"},
            ]}
            scheduler=_RollingScheduler(run,cfg)
            with patch("cleangene.workers.user_queue_snapshot",return_value=snapshot): scheduler.refresh()
            self.assertEqual(scheduler.active_indices("preprocess"),{1,2})
            self.assertEqual(scheduler.stage_queue("preprocess"),(1,1))
            output=StringIO()
            with contextlib.redirect_stdout(output): scheduler.progress("preprocess",[0,1,2],"CleanGene preprocess")
            text=output.getvalue()
            self.assertIn("running=1",text)
            self.assertIn("slurm_pending=1",text)

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
            text=output.getvalue()
            self.assertRegex(text,r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \| step=CleanGene preprocess")
            self.assertIn("total_completed=1",text)
            self.assertIn("current_step_completed=1/3",text)
            self.assertIn("total_submitted=3",text)
            self.assertIn("not_submitted_yet=0",text)

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
            output=StringIO()
            with patch("cleangene.workers._controller_pipeline") as pipeline, patch("cleangene.workers._run_single_job") as singles, contextlib.redirect_stdout(output): slurm_controller(run)
            pipeline.assert_called_once_with(run,True)
            self.assertFalse(any("resolve_groups" in str(c) for c in singles.call_args_list))
            self.assertIn("controller_started",output.getvalue())
            self.assertIn("stages=kraken_db_setup -> preprocess -> resolve_groups",output.getvalue())

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
                        if stage=="reduce": write_tsv(run/"results"/"groups"/name/"cleaned_pangenome.tsv",["Gene"],[])
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
            self.assertEqual(read_tsv(root/"cleaned_pangenome.tsv"),rows)
            self.assertTrue(_index_done(run,"reduce",0))

    def test_reduce_done_without_primary_output_is_incomplete(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["g",2,"small"]]); atomic_json(run/"state"/"reduce"/"g.done.json",{})
            self.assertFalse(_index_done(run,"reduce",0))

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
