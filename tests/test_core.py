import tempfile, unittest
from pathlib import Path
from io import StringIO
import contextlib
from cleangene.manifest import groups, load_manifest
from cleangene.fasta import assembly_metrics
from cleangene.pangenome import normalize_panaroo, present, select_rows
from cleangene.slurm import array_chunks, sbatch_cmd, submit, submit_chunked_arrays
from cleangene.util import atomic_json, write_tsv
from cleangene.workers import ensure_kraken2_db, manifest_pangenome_dir, parse_kraken_report
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
            manifest.write_text("isolate_id\torganism\tR1\tR2\n" + "\n".join(f"i{x}\tSpecies\tr1.fq\tr2.fq" for x in range(3)) + "\n")
            cfg={"SLURM_ACCOUNT":"","SLURM_PARTITION":"","SLURM_MAX_PARALLEL":"100","SLURM_ARRAY_CHUNK_SIZE":"2","SLURM_MAX_OUTSTANDING_CHUNKS":"1","SLURM_CPUS":"1","SLURM_MEM":"1G","SLURM_TIME":"1:00:00","GROUP_ORCHESTRATOR_CPUS":"1","GROUP_ORCHESTRATOR_MEM":"1G","GROUP_ORCHESTRATOR_TIME":"1:00:00","TAXONOMY_MODE":"kraken2","KRAKEN2_DB":"","KRAKEN2_DB_CPUS":"1","KRAKEN2_DB_MEM":"1G","KRAKEN2_DB_TIME":"1:00:00"}
            run=make_run(manifest,Path(d),cfg,"r")
            buf=StringIO()
            with contextlib.redirect_stdout(buf):
                slurm(run,cfg,True)
            out=buf.getvalue()
            self.assertIn("cg-kraken_db_setup",out)
            self.assertIn("--array 0-1%100",out)
            self.assertIn("--array 2%100",out)
            self.assertIn("--dependency afterok:DRYRUN",out)

    def test_preprocess_cannot_build_missing_kraken_db(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"runs"/"r"; (run/"provenance").mkdir(parents=True)
            cfg={"TAXONOMY_MODE":"kraken2","KRAKEN2_DB":"","KRAKEN2_AUTO_DOWNLOAD":"true","KRAKEN2_DATABASE_SIZE":"standard-8"}
            rows=[{"isolate_id":"i1","group_id":"g","grouping_source":"manifest_group_id","R1":"r1","R2":"r2"}]
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","grouping_source","R1","R2"],rows)
            with self.assertRaisesRegex(SystemExit,"kraken_db_setup"):
                ensure_kraken2_db(run,cfg,rows)

if __name__ == "__main__": unittest.main()
