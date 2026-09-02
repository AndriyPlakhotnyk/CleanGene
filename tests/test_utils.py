from __future__ import annotations
import contextlib, json, tempfile, unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from cleangene.cli import main
from cleangene.defaults import DEFAULTS
from cleangene.diagnostics import infer_failure_stage
from cleangene.downstream import differential_genes, get_operon, get_samples, get_variants, make_itol, matrix_path, resolve_organism
from cleangene.util import atomic_json, read_tsv, write_tsv

class CleanGeneUtilsTests(unittest.TestCase):
    def test_utils_prefer_primary_cleaned_pangenome(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); root=run/"results"/"groups"/"Species_one"
            write_tsv(root/"cleaned_pangenome.tsv",["Gene","i1"],[["primary",1]])
            write_tsv(root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv",["Gene","i1"],[["legacy",1]])
            self.assertEqual(matrix_path(run,"Species one"),root/"cleaned_pangenome.tsv")

    def test_diagnostic_failure_stage_tracks_sampling_and_assembler(self):
        self.assertEqual(infer_failure_stage(0,1,20,20,0,0,0,0,False,False)[0],"shovill_seqkit_sampling")
        self.assertEqual(infer_failure_stage(0,1,20,20,20,0,0,0,False,False)[0],"spades_graph_construction")
        self.assertEqual(infer_failure_stage(0,1,20,20,20,20,0,0,False,False)[0],"spades_contig_resolution")
        self.assertEqual(infer_failure_stage(0,1,20,20,20,20,20,0,False,False)[0],"shovill_post_assembly")
        self.assertEqual(infer_failure_stage(0,1,20,20,20,20,20,20,True,False)[0],"prokka")
    def make_run(self, root: Path) -> Path:
        run=root/"runs"/"test_run"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir(); matrix=run/"results"/"groups"/"Species_one"/"03_read_validation"/"validated_gene_presence_absence.binary.tsv"
        atomic_json(run/"provenance"/"resolved_config.json",{"SLURM_ACCOUNT":"","SLURM_PARTITION":"","UTILS_CPUS":"2","UTILS_MEM":"4G","UTILS_TIME":"01:00:00","UTILS_VARIANT_CPUS":"4","UTILS_VARIANT_MEM":"8G","UTILS_VARIANT_TIME":"02:00:00"})
        write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id"],[["i1","Species one"],["i2","Species one"],["i3","Species one"],["i4","Species one"]])
        write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],[["Species one",4,"small"]])
        write_tsv(matrix,["Gene","i1","i2","i3","i4"],[["gA",1,1,0,0],["gB",1,0,1,0],["gC",0,0,1,1]])
        return run

    def test_get_samples_writes_calls_and_lists(self):
        with tempfile.TemporaryDirectory() as d:
            run=self.make_run(Path(d)); out=Path(d)/"out"; out.mkdir()
            get_samples({"run_dir":str(run),"output_dir":str(out),"organism":"Species one","genes":["gA","gB"],"samples":[],"status":"both","match":"all"})
            calls=read_tsv(out/"gene_presence_absence.tsv"); present=read_tsv(out/"samples_present.tsv"); absent=read_tsv(out/"samples_absent.tsv")
            self.assertEqual([r["isolate_id"] for r in present],["i1"]); self.assertEqual(len(absent),3); self.assertEqual(calls[0]["gB"],"1")

    def test_differential_genes_compares_two_cohorts(self):
        with tempfile.TemporaryDirectory() as d, patch("cleangene.downstream.plot_binary_heatmap"):
            run=self.make_run(Path(d)); out=Path(d)/"out"; out.mkdir()
            differential_genes({"run_dir":str(run),"output_dir":str(out),"organism":"Species one","group_a":["i1","i2"],"group_b":["i3","i4"],"max_q_value":1,"min_prevalence_difference":0,"top":10})
            rows={r["Gene"]:r for r in read_tsv(out/"differential_genes.tsv")}
            self.assertEqual(rows["gA"]["group_a_present"],"2"); self.assertEqual(rows["gA"]["group_b_present"],"0")

    def test_operon_assigns_same_id_to_same_pattern(self):
        with tempfile.TemporaryDirectory() as d, patch("cleangene.downstream.plot_binary_heatmap"):
            run=self.make_run(Path(d)); out=Path(d)/"out"; out.mkdir()
            get_operon({"run_dir":str(run),"output_dir":str(out),"operons":[{"name":"abc","genes":["gA"],"organism":"Species one","samples":[]}]})
            rows=read_tsv(out/"operon_calls.tsv"); by_iso={r["isolate_id"]:r for r in rows}
            self.assertEqual(by_iso["i3"]["operon_id"],by_iso["i4"]["operon_id"]); self.assertEqual(by_iso["i1"]["operon_id"],by_iso["i2"]["operon_id"]); self.assertNotEqual(by_iso["i1"]["operon_id"],by_iso["i3"]["operon_id"])

    def test_itol_writes_binary_and_operon_strip(self):
        with tempfile.TemporaryDirectory() as d:
            run=self.make_run(Path(d)); out=Path(d)/"out"; out.mkdir(); operon=Path(d)/"operon.tsv"
            write_tsv(operon,["isolate_id","operon_id"],[["i1","abc_V001"],["i2","abc_V002"]])
            make_itol({"run_dir":str(run),"output_dir":str(out),"organism":"Species one","genes":["gA"],"samples":[],"operon":str(operon),"variants":"","color_scheme":"muted","custom_colors":{}})
            self.assertTrue((out/"itol_gene_presence_absence.txt").read_text().startswith("DATASET_BINARY")); self.assertTrue((out/"itol_operon_types.txt").read_text().startswith("DATASET_COLORSTRIP"))

    def test_missing_organism_lists_run_organisms(self):
        with tempfile.TemporaryDirectory() as d:
            run=self.make_run(Path(d))
            with self.assertRaisesRegex(SystemExit,"Available organisms: Species one"): resolve_organism(run,None,[])

    def test_variant_worker_reports_read_metrics_and_location_status(self):
        with tempfile.TemporaryDirectory() as d, \
             patch("cleangene.workers.retained_rows",return_value=[{"isolate_id":"i1","R1":"r1","R2":"r2","assembly":"","gff":""}]), \
             patch("cleangene.downstream._gene_references",return_value=([("CG00000001","ACGT")],[{"reference_id":"CG00000001","Gene":"gA","feature_type":"target","parent_gene":"","flank_offset":"","annotation":"gA"}])), \
             patch("cleangene.downstream.run"), patch("cleangene.evidence.map_reads"), \
             patch("cleangene.evidence.coverage",return_value={"CG00000001":{"mapped_reads":10,"covered_bases":4,"breadth":1.0}}), \
             patch("cleangene.evidence.align_identity",return_value={"CG00000001":{"identity":1.0,"identical_positions":4}}), \
             patch("cleangene.downstream._variant_counts",return_value={"CG00000001":{"snps":1,"mnps":0,"insertions":0,"deletions":0,"inserted_bases":0,"deleted_bases":0}}), \
             patch("cleangene.downstream._plot_variant_alignment"):
            run_dir=self.make_run(Path(d)); out=Path(d)/"out"; out.mkdir()
            def fake_consensus(reference,bam,prefix,*args):
                path=prefix.with_suffix(".consensus.fasta"); path.write_text(">CG00000001\nACGT\n"); return path
            with patch("cleangene.evidence.consensus",side_effect=fake_consensus):
                get_variants({"run_dir":str(run_dir),"output_dir":str(out),"organism":"Species one","genes":["gA"],"samples":["i1"],"cpus":1,"min_similarity":95,"flanking_genes":0})
            row=read_tsv(out/"gene_variants.tsv")[0]
            self.assertEqual(row["variant_id"],"gA_V001"); self.assertEqual(row["percent_identity"],"100.0"); self.assertEqual(row["snps"],"1"); self.assertEqual(row["location_status"],"no_assembly")

    def test_utils_cli_prints_submission_messages_and_uses_sbatch(self):
        with tempfile.TemporaryDirectory() as d, patch("cleangene.utils_cli.submit",return_value="123") as submit_job:
            run=self.make_run(Path(d)); stdout=StringIO()
            with contextlib.redirect_stdout(stdout):
                code=main(["utils","get-samples","--run-dir",str(run),"--analysis-name","query","--organism","Species one","--genes","gA"])
            self.assertEqual(code,0); text=stdout.getvalue(); self.assertIn("Welcome to CleanGene Utils!",text); self.assertIn(f"Located run: {run}",text); self.assertIn("Getting ready to submit",text); self.assertIn("Analysis submitted. Please find logs in",text)
            command=submit_job.call_args.args[0]; self.assertEqual(command[0],"sbatch"); self.assertIn("_utils_worker",command[-1])
            request=json.loads((run/"results"/"utils"/"query"/"request.json").read_text()); self.assertEqual(request["slurm_job_id"],"123")

    def test_checkm2_utility_submits_setup_array_and_merge_resources(self):
        with tempfile.TemporaryDirectory() as d, patch("cleangene.utils_cli.submit",side_effect=["setup","array","merge"]) as submit_job:
            run=self.make_run(Path(d)); sample=run/"results"/"sample_data"/"i1"; assembly=sample/"assembly"/"contigs.fasta"; assembly.parent.mkdir(parents=True); assembly.write_text(">c\nACGT\n")
            write_tsv(run/"results"/"cohort"/"isolate_qc.tsv",["isolate_id","group_id","assembly","PASS/FAIL","Notes","excluded","reason","top_species","contamination_pct","trimmed_read_length","mean_base_quality","sequencing_coverage","contigs","n50","gff","checkm2_completeness","checkm2_contamination"],[["i1","Species one",assembly,"PASS","All evaluated QC criteria passed","0","","","",150,35,40,80,75000,"gff","",""]])
            cfg=dict(DEFAULTS); cfg.update({"CHECKM2_MODE":"off","CHECKM2_CPUS":"16","CHECKM2_MEM":"128G","CHECKM2_TIME":"24:00:00","CHECKM2_POSTHOC_CPUS":"1","CHECKM2_POSTHOC_MEM":"32G","CHECKM2_POSTHOC_TIME":"02:00:00","CHECKM2_POSTHOC_MAX_INFLIGHT":"8"})
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            with contextlib.redirect_stdout(StringIO()):
                self.assertEqual(main(["utils","checkm2","--run-dir",str(run),"--analysis-name","posthoc","--samples","i1"]),0)
            commands=[call.args[0] for call in submit_job.call_args_list]
            self.assertIn("cg-checkm2_posthoc_setup",commands[0])
            self.assertEqual(commands[0][commands[0].index("--cpus-per-task")+1],"16")
            self.assertEqual(commands[0][commands[0].index("--mem")+1],"128G")
            self.assertIn("cg-util-checkm2",commands[1])
            self.assertEqual(commands[1][commands[1].index("--cpus-per-task")+1],"1")
            self.assertEqual(commands[1][commands[1].index("--mem")+1],"32G")
            self.assertIn("--dependency",commands[1])
            self.assertIn("cg-util-checkm2-merge",commands[2])
            request=json.loads((run/"results"/"utils"/"posthoc"/"checkm2"/"request.json").read_text())
            self.assertEqual(request["setup_slurm_job_id"],"setup")
            self.assertEqual(request["array_slurm_job_id"],"array")
            self.assertEqual(request["merge_slurm_job_id"],"merge")

    def test_checkm2_posthoc_merge_enriches_cohort_without_rewriting_original_qc(self):
        from cleangene.downstream import checkm2_posthoc_merge
        with tempfile.TemporaryDirectory() as d:
            run=self.make_run(Path(d)); out=run/"results"/"utils"/"posthoc"/"checkm2"; out.mkdir(parents=True)
            write_tsv(out/"tasks.tsv",["isolate_id","group_id","assembly"],[["i1","Species one","asm.fa"]])
            iso=out/"isolates"/"i1"; iso.mkdir(parents=True)
            write_tsv(iso/"checkm2_qc.tsv",["isolate_id","group_id","assembly","checkm2_completeness","checkm2_contamination","checkm2_status","checkm2_notes"],[["i1","Species one","asm.fa","74","1","FAIL","FAIL: completeness=74% is below fail minimum 80%"]])
            cohort=run/"results"/"cohort"; write_tsv(cohort/"isolate_qc.tsv",["isolate_id","group_id","assembly","PASS/FAIL","Notes","excluded","reason","top_species","contamination_pct","trimmed_read_length","mean_base_quality","sequencing_coverage","contigs","n50","gff","checkm2_completeness","checkm2_contamination"],[["i1","Species one","asm.fa","PASS","All evaluated QC criteria passed","0","","","",150,35,40,80,75000,"i1.gff","",""]])
            checkm2_posthoc_merge({"run_dir":str(run),"output_dir":str(out)})
            row=read_tsv(cohort/"isolate_qc.tsv")[0]
            self.assertEqual(row["PASS/FAIL"],"PASS")
            self.assertEqual(row["excluded"],"0")
            self.assertEqual(row["checkm2_completeness"],"74")
            self.assertEqual(row["checkm2_posthoc_status"],"FAIL")
            self.assertEqual(row["PASS/FAIL_with_checkm2"],"FAIL")
            self.assertIn("Panaroo membership was not changed",(out/"summary.txt").read_text())

    def test_diagnostic_cli_submits_with_dedicated_resources(self):
        with tempfile.TemporaryDirectory() as d, patch("cleangene.utils_cli.submit",return_value="456") as submit_job:
            run=self.make_run(Path(d))
            with contextlib.redirect_stdout(StringIO()): main(["utils","diagnose-call","--run-dir",str(run),"--analysis-name","diag","--organism","Species one","--genes","gA","--samples","i1"])
            command=submit_job.call_args.args[0]; self.assertEqual(command[command.index("--cpus-per-task")+1],"16"); self.assertEqual(command[command.index("--mem")+1],"128G")
            request=json.loads((run/"results"/"utils"/"diag"/"request.json").read_text()); self.assertEqual(request["utility"],"diagnose_call"); self.assertEqual(request["samples"],["i1"])

    def test_core_run_prints_submission_messages(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); manifest=root/"manifest.tsv"; manifest.write_text("isolate_id\torganism\tR1\tR2\ni1\tSpecies one\t/r1.fastq.gz\t/r2.fastq.gz\n")
            config=root/"off.env"; config.write_text("CHECKM2_MODE=off\n")
            stdout=StringIO()
            with contextlib.redirect_stdout(stdout): main(["run","--manifest",str(manifest),"--analysis-root",str(root),"--run-id","messages","--dry-run","--config",str(config)])
            text=stdout.getvalue(); self.assertTrue(text.startswith("+")); self.assertIn("Cleanse thy pangenome, my liege!",text); self.assertIn(f"Run directory: {root/'runs'/'messages'}",text); self.assertIn("Getting ready to submit",text); self.assertIn("CleanGene run created: messages",text); self.assertIn("Controller submitted: DRYRUN",text); self.assertIn("Welcome to CleanGene, Your Grace.",text); self.assertLess(text.index("Cleanse thy pangenome, my liege!"),text.index("Welcome to CleanGene, Your Grace.")); self.assertLess(text.index("Welcome to CleanGene, Your Grace."),text.index("Run directory:"))

    def test_resume_submits_without_login_side_legacy_scans(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"runs"/"resume"; (run/"provenance").mkdir(parents=True); (run/"state").mkdir()
            cfg=dict(DEFAULTS); cfg["CHECKM2_MODE"]="off"; atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","organism"],[["i1","Species one","Species one"]])
            write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],[["Species one","i1"]])
            stdout=StringIO()
            with patch("cleangene.cli.slurm",return_value="123") as submit_job, patch("cleangene.cli.invalidate_legacy_identity_metrics",side_effect=AssertionError("slow identity scan ran")), patch("cleangene.cli.invalidate_legacy_isolate_qc",side_effect=AssertionError("slow QC scan ran")), contextlib.redirect_stdout(stdout):
                self.assertEqual(main(["resume","--run-dir",str(run)]),0)
            self.assertEqual(submit_job.call_count,1); text=stdout.getvalue(); self.assertTrue(text.startswith("+")); self.assertIn("Cleanse thy pangenome, my liege!",text); self.assertIn("legacy checks will run inside the controller job",text); self.assertIn("CleanGene run created: resume",text); self.assertIn("Controller submitted: 123",text); self.assertIn("Welcome to CleanGene, Your Grace.",text); self.assertLess(text.index("Welcome to CleanGene, Your Grace."),text.index("Run directory:"))

if __name__=="__main__": unittest.main()
