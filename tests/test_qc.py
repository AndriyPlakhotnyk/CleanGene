import contextlib, json, tempfile, unittest
from io import StringIO
from pathlib import Path

from cleangene.cli import check, invalidate_legacy_isolate_qc
from cleangene.qc import (QC_OUTPUT_FIELDS, THRESHOLD_DEFAULTS, classify_isolate_qc,
    parse_checkm2_report, prepare_qc_provenance, read_metrics, resolve_threshold_rows,
    validate_thresholds)
from cleangene.util import atomic_json, read_tsv, write_tsv
from cleangene.workers import retained_rows
from unittest.mock import patch


class IsolateQCTests(unittest.TestCase):
    def setUp(self):
        self.thresholds=validate_thresholds(THRESHOLD_DEFAULTS)

    def classify(self,**updates):
        values=dict(thresholds=self.thresholds,expected_organism="Escherichia coli",top_species="  escherichia   COLI ",
            kraken_contamination=5.0,read_length=120.0,mean_quality=30.0,coverage=20.0,
            contigs=300.0,n50=25000.0,completeness=90.0,checkm2_contamination=5.0,
            checkm2_mode="required",internal_pangenome=True,external_pangenome=False,
            assembly_present=True,gff_present=True)
        values.update(updates); return classify_isolate_qc(**values)

    def test_all_pass_boundaries_are_pass(self):
        result=self.classify()
        self.assertEqual(result["PASS/FAIL"],"PASS")
        self.assertEqual(result["Notes"],"All evaluated QC criteria passed")
        self.assertEqual(result["excluded"],"0")

    def test_warning_boundaries_belong_to_warning(self):
        cases={"completeness":80.0,"checkm2_contamination":10.0,"contigs":1000.0,"n50":5000.0,"coverage":10.0}
        for key,value in cases.items():
            with self.subTest(metric=key):
                result=self.classify(**{key:value}); self.assertEqual(result["PASS/FAIL"],"WARNING"); self.assertEqual(result["excluded"],"0")

    def test_values_beyond_fail_boundaries_fail(self):
        cases={"kraken_contamination":5.0001,"completeness":79.999,"checkm2_contamination":10.001,"contigs":1001.0,"n50":4999.0,"coverage":9.999}
        for key,value in cases.items():
            with self.subTest(metric=key):
                result=self.classify(**{key:value}); self.assertEqual(result["PASS/FAIL"],"FAIL"); self.assertEqual(result["excluded"],"1")

    def test_read_length_and_quality_warn_without_default_fail(self):
        result=self.classify(read_length=1.0,mean_quality=1.0)
        self.assertEqual(result["PASS/FAIL"],"WARNING")
        self.assertIn("read_length=1 bp is below pass minimum 120 bp",result["Notes"])
        self.assertIn("mean_base_quality=1 is below pass minimum 30",result["Notes"])

    def test_numeric_read_fail_overrides_create_normal_bands(self):
        thresholds=dict(self.thresholds); thresholds["qc_min_read_length_fail"]=100; thresholds["qc_min_mean_base_quality_fail"]=20
        self.assertEqual(self.classify(thresholds=thresholds,read_length=99)["PASS/FAIL"],"FAIL")
        self.assertEqual(self.classify(thresholds=thresholds,mean_quality=19)["PASS/FAIL"],"FAIL")

    def test_taxonomy_unavailable_warns_mismatch_fails_and_absent_expected_is_ignored(self):
        self.assertEqual(self.classify(top_species="")["PASS/FAIL"],"WARNING")
        self.assertEqual(self.classify(top_species="Klebsiella pneumoniae")["PASS/FAIL"],"FAIL")
        self.assertEqual(self.classify(expected_organism="",top_species="")["PASS/FAIL"],"PASS")

    def test_checkm2_off_is_informational_not_a_warning(self):
        result=self.classify(checkm2_mode="off",completeness=None,checkm2_contamination=None)
        self.assertEqual(result["PASS/FAIL"],"PASS"); self.assertIn("INFO: CheckM2 was not evaluated",result["Notes"])

    def test_internal_missing_assembly_or_gff_fails(self):
        self.assertEqual(self.classify(assembly_present=False,gff_present=None)["PASS/FAIL"],"FAIL")
        self.assertEqual(self.classify(gff_present=False)["PASS/FAIL"],"FAIL")

    def test_external_pangenome_without_assembly_warns(self):
        result=self.classify(expected_organism="",top_species="",coverage=None,contigs=None,n50=None,
            completeness=None,checkm2_contamination=None,internal_pangenome=False,external_pangenome=True,
            assembly_present=False,gff_present=None)
        self.assertEqual(result["PASS/FAIL"],"WARNING"); self.assertEqual(result["excluded"],"0")
        self.assertIn("external pangenome was supplied without an assembly",result["Notes"])

    def test_multiple_notes_are_deterministic_and_include_values(self):
        result=self.classify(coverage=8.7,contigs=450)
        self.assertEqual(result["PASS/FAIL"],"FAIL")
        self.assertIn("FAIL: coverage=8.7x is below fail minimum 10x",result["Notes"])
        self.assertIn("WARNING: contigs=450 exceeds pass maximum 300",result["Notes"])
        self.assertLess(result["Notes"].index("coverage=8.7x"),result["Notes"].index("contigs=450"))

    def test_user_exclusion_is_fail(self):
        result=self.classify(user_exclusion=True)
        self.assertEqual(result["PASS/FAIL"],"FAIL"); self.assertEqual(result["reason"],"user_excluded")

    def test_profile_precedence_and_provenance_copy(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); profile=root/"profiles.tsv"
            profile.write_text("scope_type\tscope_value\tqc_max_contigs_pass\norganism\tSpecies A\t250\ngroup_id\tg1\t200\n")
            rows=[{"isolate_id":"i1","organism":"Species A","group_id":"g1","qc_max_contigs_pass":"150"},
                  {"isolate_id":"i2","organism":"Species A","group_id":"g1"},
                  {"isolate_id":"i3","organism":" species   a ","group_id":"g2"}]
            resolved=resolve_threshold_rows(rows,{"QC_PROFILE_FILE":str(profile)},profile)
            self.assertEqual([r["qc_max_contigs_pass"] for r in resolved],["150","200","250"])
            self.assertEqual([r["qc_profile_source"] for r in resolved],["manifest:i1","group_id:g1","organism: species   a "])
            run=root/"run"; cfg=prepare_qc_provenance(run,rows,{"QC_PROFILE_FILE":str(profile)})
            self.assertEqual(Path(cfg["QC_PROFILE_FILE"]),run/"provenance"/"qc_profile.tsv")
            self.assertTrue((run/"provenance"/"qc_thresholds.tsv").is_file())

    def test_blank_manifest_override_inherits(self):
        resolved=resolve_threshold_rows([{"isolate_id":"i","group_id":"g","qc_max_contigs_pass":""}],{"QC_MAX_CONTIGS_PASS":"321"})[0]
        self.assertEqual(resolved["qc_max_contigs_pass"],"321")

    def test_invalid_profiles_numbers_negatives_and_ordering_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); row=[{"isolate_id":"i","group_id":"g"}]
            for content in (
                "scope_type\tscope_value\norganism\tA\norganism\t a \n",
                "scope_type\tscope_value\tqc_max_contigs_pass\ngroup_id\tg\tabc\n",
                "scope_type\tscope_value\tqc_max_contigs_pass\ngroup_id\tg\t-1\n",
                "scope_type\tscope_value\tqc_max_contigs_pass\ngroup_id\tg\tnan\n",
                "scope_type\tscope_value\tqc_typo\ngroup_id\tg\t5\n",
            ):
                profile=root/"p.tsv"; profile.write_text(content)
                with self.subTest(content=content), self.assertRaises(SystemExit): resolve_threshold_rows(row,{},profile)
        with self.assertRaises(SystemExit): resolve_threshold_rows([{"isolate_id":"i","group_id":"g","qc_max_contigs_pass":"1001"}],{})

    def test_fastq_metrics_use_weighted_phred_and_valid_fastp_lengths(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); r1=root/"r1.fastq.gz"; r2=root/"r2.fastq.gz"
            r1.write_text("@a\nAAAA\n+\nIIII\n@b\nCCCC\n+\nIIII\n"); r2.write_text("@c\nGG\n+\n!!\n")
            fastp=root/"fastp.json"; fastp.write_text(json.dumps({"summary":{"after_filtering":{"total_bases":12,"read1_mean_length":6,"read2_mean_length":3}}}))
            metrics=read_metrics(r1,r2,fastp)
            self.assertEqual(metrics["trimmed_read_length"],3); self.assertEqual(metrics["total_bases"],12); self.assertEqual(metrics["mean_base_quality"],32)

    def test_checkm2_report_parser(self):
        with tempfile.TemporaryDirectory() as d:
            report=Path(d)/"quality_report.tsv"; report.write_text("Name\tCompleteness\tContamination\ni\t90.5\t4.25\n")
            self.assertEqual(parse_checkm2_report(report),(90.5,4.25))

    def test_warning_is_retained_and_fail_is_excluded(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); atomic_json(run/"provenance"/"resolved_config.json",{})
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id"],[["warn","g"],["fail","g"]])
            for isolate,status,excluded in (("warn","WARNING",0),("fail","FAIL",1)):
                qc=run/"results"/"groups"/"g"/"01_isolates"/isolate/"qc.tsv"
                write_tsv(qc,["isolate_id","group_id","excluded","PASS/FAIL","assembly","gff","R1","R2"],[[isolate,"g",excluded,status,"","","r1","r2"]])
            self.assertEqual([row["isolate_id"] for row in retained_rows(run,"g")],["warn"])

    def test_legacy_qc_marker_is_invalidated(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); qc=run/"results"/"groups"/"g"/"01_isolates"/"i"/"qc.tsv"
            write_tsv(qc,["isolate_id","excluded"],[['i',0]]); marker=run/"state"/"preprocess"/"i.done.json"; atomic_json(marker,{})
            self.assertEqual(invalidate_legacy_isolate_qc(run),1); self.assertFalse(marker.exists())

    def test_required_output_headers_are_exact(self):
        self.assertEqual(QC_OUTPUT_FIELDS,("PASS/FAIL","Notes","trimmed_read_length","mean_base_quality","sequencing_coverage","checkm2_completeness","checkm2_contamination","qc_profile_source"))

    def test_check_requires_checkm2_executable_and_database(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); manifest=root/"manifest.tsv"; manifest.write_text("isolate_id\tgroup_id\tR1\tR2\ni\tg\tr1\tr2\n")
            config=root/"required.env"; config.write_text("ASSEMBLER=off\nTAXONOMY_MODE=off\nCHECKM2_MODE=required\nCHECKM2_DB=\n")
            args=type("Args",(),{"manifest":manifest,"config":config})()
            with patch("cleangene.cli.command_exists",return_value=True),patch("cleangene.cli.resolve_checkm2_executable",return_value=root/"checkm2"),contextlib.redirect_stdout(StringIO()): self.assertEqual(check(args),0)
            database=root/"uniref100.KO.1.dmnd"; database.write_text("db"); config.write_text(f"ASSEMBLER=off\nTAXONOMY_MODE=off\nCHECKM2_MODE=required\nCHECKM2_DB={database}\n")
            with patch("cleangene.cli.command_exists",return_value=True),patch("cleangene.cli.resolve_checkm2_executable",return_value=root/"checkm2"),contextlib.redirect_stdout(StringIO()): self.assertEqual(check(args),0)
            config.write_text("ASSEMBLER=off\nTAXONOMY_MODE=off\nCHECKM2_MODE=off\nCHECKM2_DB=\n")
            with patch("cleangene.cli.command_exists",return_value=False),contextlib.redirect_stdout(StringIO()): self.assertEqual(check(args),0)


if __name__=="__main__": unittest.main()
