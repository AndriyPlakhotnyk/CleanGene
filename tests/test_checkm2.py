import multiprocessing
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from cleangene.checkm2 import CheckM2DbError, CheckM2DbNotReady, EXPECTED_CHECKM2_DB_NAME, CHECKM2_COMMAND_SCHEMA_VERSION, bundled_test_genome, checkm2_database_root, checkm2_named_input_link, checkm2_predict_capabilities, checkm2_predict_command, checkm2_runtime_marker, parse_checkm2_quality_report, resolve_checkm2_db, validate_checkm2_db
from cleangene.defaults import DEFAULTS
from cleangene.tools import ToolResolutionError, resolve_checkm2_executable
from cleangene.cli import make_run
from cleangene.util import atomic_json, load_json, read_tsv, write_tsv
from cleangene.workers import checkm2_db_setup, preprocess, slurm_controller


def _write_db(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("diamond db\n")

def _write_executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("""#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("CheckM2 version 1.1.0")
elif args[:2] == ["predict", "--help"] or args == ["predict", "-h"]:
    print("--input --output-directory --database_path --threads --remove_intermediates --force")
elif args and args[0] == "testrun":
    sys.exit(0)
elif args and args[0] == "predict":
    out = pathlib.Path(args[args.index("--output-directory") + 1])
    input_path = pathlib.Path(args[args.index("--input") + 1])
    out.mkdir(parents=True, exist_ok=True)
    (out / "quality_report.tsv").write_text(f"Name\\tCompleteness\\tContamination\\n{input_path.name.split('.')[0]}\\t99\\t1\\n")
else:
    print("unsupported", args, file=sys.stderr)
    sys.exit(2)
""")
    path.chmod(0o755)
    return path


def _concurrent_checkm2_worker(root: str, calls_dir: str, results_dir: str, errors_dir: str) -> None:
    cfg = {"CHECKM2_DATABASE_ROOT": root, "CHECKM2_AUTO_DOWNLOAD": "true", "CHECKM2_EXECUTABLE": str(Path(root).parent / "bin" / "checkm2")}

    def runner(command):
        Path(calls_dir, f"call-{time.time_ns()}").write_text("\t".join(command))
        time.sleep(0.2)
        _write_db(Path(command[command.index("--path") + 1]) / "CheckM2_database" / EXPECTED_CHECKM2_DB_NAME)

    try:
        result = resolve_checkm2_db(cfg, allow_download=True, runner=runner)
        Path(results_dir, f"result-{time.time_ns()}").write_text(str(result.path))
    except Exception as error:
        Path(errors_dir, f"error-{time.time_ns()}").write_text(str(error))


class CheckM2DatabaseTests(unittest.TestCase):
    def test_existing_shared_db_is_reused_and_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            exe=_write_executable(Path(d)/"bin"/"checkm2")
            run = Path(d) / "run"; cfg = {"CHECKM2_MODE": "required", "CHECKM2_DATABASE_ROOT": str(Path(d) / "checkm2"), "CHECKM2_AUTO_DOWNLOAD": "true", "CHECKM2_EXECUTABLE": str(exe)}
            db = checkm2_database_root(cfg) / "CheckM2_database" / EXPECTED_CHECKM2_DB_NAME; _write_db(db)
            atomic_json(run / "provenance" / "resolved_config.json", cfg)
            write_tsv(run / "provenance" / "manifest.tsv", ["isolate_id", "group_id", "R1", "R2"], [["i1", "g", "r1", "r2"]])
            with patch("cleangene.workers.bundled_test_genome", return_value=db):
                checkm2_db_setup(run)
            resolved = load_json(run / "provenance" / "resolved_config.json")
            marker = load_json(run / "state" / "checkm2_db_setup.done.json")
            self.assertEqual(resolved["CHECKM2_DB"], str(db.resolve()))
            self.assertEqual(resolved["CHECKM2_EXECUTABLE"], str(exe.resolve()))
            self.assertEqual(resolved["CHECKM2_VERSION"], "CheckM2 version 1.1.0")
            self.assertEqual(marker["source"], "shared_existing")
            self.assertTrue(marker["runtime_verified"])
            self.assertTrue(checkm2_runtime_marker(cfg).is_file())

    def test_pinned_cli_uses_underscore_cleanup_option(self):
        with tempfile.TemporaryDirectory() as d:
            exe=_write_executable(Path(d)/"bin"/"checkm2")
            caps=checkm2_predict_capabilities(exe)
            command=checkm2_predict_command(exe,Path(d)/"i.fna",Path(d)/"out",Path(d)/EXPECTED_CHECKM2_DB_NAME,1,caps)
            self.assertIn("--remove_intermediates",command)
            self.assertNotIn("--remove-intermediates",command)

    def test_runtime_marker_invalidates_when_schema_changes(self):
        with tempfile.TemporaryDirectory() as d:
            exe=_write_executable(Path(d)/"bin"/"checkm2")
            cfg={"CHECKM2_DATABASE_ROOT":str(Path(d)/"checkm2"),"CHECKM2_EXECUTABLE":str(exe)}
            db=checkm2_database_root(cfg)/"CheckM2_database"/EXPECTED_CHECKM2_DB_NAME; _write_db(db)
            from cleangene.checkm2 import record_checkm2_runtime_verified, checkm2_runtime_is_verified
            record_checkm2_runtime_verified(cfg,db,exe,"CheckM2 version 1.1.0")
            marker=checkm2_runtime_marker(cfg)
            data=load_json(marker); data["schema"]=CHECKM2_COMMAND_SCHEMA_VERSION-1; atomic_json(marker,data)
            self.assertFalse(checkm2_runtime_is_verified(cfg,db,exe,"CheckM2 version 1.1.0"))

    def test_runtime_smoke_test_failure_stops_before_preprocess(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); exe=_write_executable(root/"bin"/"checkm2"); run=root/"run"
            cfg={"CHECKM2_MODE":"required","CHECKM2_DATABASE_ROOT":str(root/"checkm2"),"CHECKM2_AUTO_DOWNLOAD":"true","CHECKM2_EXECUTABLE":str(exe)}
            db=checkm2_database_root(cfg)/"CheckM2_database"/EXPECTED_CHECKM2_DB_NAME; _write_db(db)
            atomic_json(run/"provenance"/"resolved_config.json",cfg)
            write_tsv(run/"provenance"/"manifest.tsv",["isolate_id","group_id","R1","R2"],[["i1","g","r1","r2"]])

            def fail_testrun(command, **kwargs):
                kwargs["stderr"].parent.mkdir(parents=True,exist_ok=True)
                kwargs["stderr"].write_text("database checksum is incompatible\n")
                raise __import__("subprocess").CalledProcessError(1,command)

            with patch("cleangene.workers.run",side_effect=fail_testrun):
                with self.assertRaisesRegex(SystemExit,"database checksum is incompatible"):
                    checkm2_db_setup(run)
            self.assertFalse((run/"state"/"checkm2_db_setup.done.json").exists())
            self.assertFalse(checkm2_runtime_marker(cfg).exists())

    def test_missing_db_auto_downloads_once(self):
        with tempfile.TemporaryDirectory() as d:
            exe=_write_executable(Path(d)/"bin"/"checkm2")
            cfg = {"CHECKM2_DATABASE_ROOT": str(Path(d) / "checkm2"), "CHECKM2_AUTO_DOWNLOAD": "true", "CHECKM2_EXECUTABLE": str(exe)}
            calls = []

            def runner(command):
                calls.append(command)
                _write_db(Path(command[command.index("--path") + 1]) / "CheckM2_database" / EXPECTED_CHECKM2_DB_NAME)

            result = resolve_checkm2_db(cfg, allow_download=True, runner=runner)
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0][0], str(exe.resolve()))
            self.assertIn("--no_write_json_db", calls[0])
            self.assertEqual(result.source, "auto_download")
            validate_checkm2_db(result.path)

    def test_missing_db_auto_download_false_fails(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = {"CHECKM2_DATABASE_ROOT": str(Path(d) / "checkm2"), "CHECKM2_AUTO_DOWNLOAD": "false"}
            with self.assertRaisesRegex(CheckM2DbNotReady, "CHECKM2_AUTO_DOWNLOAD=false"):
                resolve_checkm2_db(cfg)

    def test_explicit_db_is_exact_and_invalid_does_not_download(self):
        with tempfile.TemporaryDirectory() as d:
            db = Path(d) / EXPECTED_CHECKM2_DB_NAME; _write_db(db)
            result = resolve_checkm2_db({"CHECKM2_DB": str(db)}, allow_download=True, runner=lambda command: (_ for _ in ()).throw(AssertionError("download should not run")))
            self.assertEqual(result.path, db.resolve())
            with self.assertRaisesRegex(CheckM2DbError, "Explicit CHECKM2_DB is invalid"):
                resolve_checkm2_db({"CHECKM2_DB": str(Path(d) / "missing.dmnd")}, allow_download=True, runner=lambda command: (_ for _ in ()).throw(AssertionError("download should not run")))

    def test_concurrent_missing_db_downloads_once(self):
        with tempfile.TemporaryDirectory() as d:
            _write_executable(Path(d)/"bin"/"checkm2")
            calls = Path(d) / "calls"; results = Path(d) / "results"; errors = Path(d) / "errors"
            calls.mkdir(); results.mkdir(); errors.mkdir()
            processes = [multiprocessing.Process(target=_concurrent_checkm2_worker, args=(str(Path(d) / "checkm2"), str(calls), str(results), str(errors))) for _ in range(2)]
            for process in processes: process.start()
            for process in processes: process.join(10)
            for process in processes: self.assertEqual(process.exitcode, 0)
            self.assertEqual(list(errors.iterdir()), [])
            self.assertEqual(len(list(calls.iterdir())), 1)
            paths = [p.read_text() for p in results.iterdir()]
            self.assertEqual(len(set(paths)), 1)
            validate_checkm2_db(Path(paths[0]))

    def test_setup_failure_stops_controller_before_preprocess(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d) / "run"; (run / "provenance").mkdir(parents=True); (run / "state").mkdir()
            cfg = {"CHECKM2_MODE": "required", "CHECKM2_EXECUTABLE": str(Path(d)/"missing"/"checkm2"), "CHECKM2_AUTO_DOWNLOAD": "false", "TAXONOMY_MODE": "off", "SLURM_USER_JOB_LIMIT": "2000", "SLURM_JOB_HEADROOM": "10", "SLURM_POLL_SECONDS": "0", "SLURM_MAX_PARALLEL": "100", "SLURM_ARRAY_CHUNK_SIZE": "500", "SLURM_ACCOUNT": "", "SLURM_PARTITION": "", "CHECKM2_CPUS": "1", "CHECKM2_MEM": "1G", "CHECKM2_TIME": "1:00:00"}
            atomic_json(run / "provenance" / "resolved_config.json", cfg)
            write_tsv(run / "provenance" / "manifest.tsv", ["isolate_id", "group_id", "R1", "R2"], [["i1", "g", "r1", "r2"]])
            write_tsv(run / "state" / "isolate_tasks.tsv", ["group_id", "isolate_id"], [["g", "i1"]])
            write_tsv(run / "state" / "group_tasks.tsv", ["group_id", "n_isolates", "group_size_class"], [["g", 1, "small"]])
            with patch("cleangene.workers._run_single_job") as single:
                with self.assertRaisesRegex(SystemExit, "Explicit checkm2 executable is invalid"):
                    slurm_controller(run)
            self.assertEqual(single.call_count,0)
            self.assertFalse((run / "results" / "sample_data" / "i1" / "qc.tsv").exists())


class CheckM2PreprocessTests(unittest.TestCase):
    def test_direct_spades_required_checkm2_then_prokka_records_sample_data(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); r1 = root / "r1.fq"; r2 = root / "r2.fq"; db = root / EXPECTED_CHECKM2_DB_NAME; exe=_write_executable(root/"bin"/"checkm2")
            r1.write_text("@r\n" + "A"*120 + "\n+\n" + "I"*120 + "\n"); r2.write_text(r1.read_text()); _write_db(db)
            manifest = root / "manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            cfg = {"TAXONOMY_MODE": "off", "ASSEMBLER": "spades", "SKIP_TRIM": "true", "CHECKM2_MODE": "required", "CHECKM2_EXECUTABLE": str(exe), "CHECKM2_DB": str(db), "READ_TRIMMING_MODE": "off", "QC_MIN_N50_PASS": "0", "QC_MIN_N50_FAIL": "0", "QC_MIN_COVERAGE_PASS": "0", "QC_MIN_COVERAGE_FAIL": "0", "PREPROCESS_USE_NODE_LOCAL_SCRATCH": "false", "COMPRESS_ASSEMBLY_OUTPUTS": "intermediates", "COMPRESS_ANNOTATION_OUTPUTS": "nonessential"}
            run = make_run(manifest, root, cfg, "r"); commands = []

            def fake_run(command, **kwargs):
                commands.append(command)
                if command[0] == "spades.py":
                    out = Path(command[command.index("-o") + 1]); out.mkdir(parents=True, exist_ok=True)
                    (out / "contigs.fasta").write_text(">c\n" + "A"*120 + "\n")
                    (out / "spades.fasta").write_text(">s\nACGT\n")
                    (out / "spades.gfa").write_text("H\tVN:Z:1.0\n")
                elif Path(command[0]).name == "checkm2":
                    out = Path(command[command.index("--output-directory") + 1]); out.mkdir(parents=True, exist_ok=True)
                    write_tsv(out / "quality_report.tsv", ["Name", "Completeness", "Contamination"], [["iso1", 99, 1]])
                elif command[0] == "prokka":
                    out = Path(command[command.index("--outdir") + 1]); prefix = command[command.index("--prefix") + 1]; out.mkdir(parents=True, exist_ok=True)
                    (out / f"{prefix}.gff").write_text("##gff-version 3\n")
                    (out / f"{prefix}.gbk").write_text("gbk\n")

            with patch("cleangene.workers.command_exists", return_value=True), patch("cleangene.workers.run", side_effect=fake_run) as runner:
                preprocess(run, 0)
            sample = run / "results" / "sample_data" / "iso1"
            self.assertFalse(any(c[0] == "fastp" for c in commands))
            self.assertFalse(any(c[0] == "shovill" for c in commands))
            self.assertIn(str(exe.resolve()), [c[0] for c in commands])
            checkm2_call=next(call for call in runner.call_args_list if Path(call.args[0][0]).name=="checkm2")
            self.assertEqual(checkm2_call.args[0][checkm2_call.args[0].index("--threads")+1],"1")
            self.assertEqual(checkm2_call.args[0][checkm2_call.args[0].index("--input")+1],str(sample/"checkm2"/"input"/"iso1.fasta"))
            self.assertIn("--remove_intermediates",checkm2_call.args[0])
            self.assertNotIn("--remove-intermediates",checkm2_call.args[0])
            self.assertEqual(checkm2_call.kwargs["env"]["TF_NUM_INTEROP_THREADS"],"1")
            self.assertEqual(checkm2_call.kwargs["env"]["OPENBLAS_NUM_THREADS"],"1")
            spades = next(c for c in commands if c[0] == "spades.py")
            self.assertIn("--only-assembler", spades)
            prokka = next(c for c in commands if c[0] == "prokka")
            self.assertEqual(Path(prokka[-1]).name, "contigs.fasta")
            qc = read_tsv(sample / "qc.tsv")[0]
            self.assertEqual(qc["PASS/FAIL"], "PASS")
            self.assertEqual(Path(qc["assembly"]), sample / "assembly" / "contigs.fasta")
            self.assertEqual(Path(qc["gff"]), sample / "annotation" / "iso1.gff")
            self.assertTrue((sample / "assembly" / "contigs.fasta").is_file())
            self.assertTrue((sample / "assembly" / "spades.fasta.gz").is_file())
            self.assertTrue((sample / "annotation" / "iso1.gff").is_file())
            self.assertTrue((sample / "state").exists() is False)
            self.assertTrue((run / "state" / "preprocess" / "iso1.done.json").is_file())

    def test_checkm2_crash_does_not_publish_fake_prokka_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); r1 = root / "r1.fq"; r2 = root / "r2.fq"; db = root / EXPECTED_CHECKM2_DB_NAME; exe=_write_executable(root/"bin"/"checkm2")
            r1.write_text("@r\n" + "A"*120 + "\n+\n" + "I"*120 + "\n"); r2.write_text(r1.read_text()); _write_db(db)
            manifest = root / "manifest.tsv"; manifest.write_text(f"isolate_id\tgroup_id\tR1\tR2\niso1\tg\t{r1}\t{r2}\n")
            cfg = {"TAXONOMY_MODE": "off", "CHECKM2_MODE": "required", "CHECKM2_EXECUTABLE": str(exe), "CHECKM2_DB": str(db), "READ_TRIMMING_MODE": "off", "QC_MIN_N50_PASS": "0", "QC_MIN_N50_FAIL": "0", "QC_MIN_COVERAGE_PASS": "0", "QC_MIN_COVERAGE_FAIL": "0", "PREPROCESS_USE_NODE_LOCAL_SCRATCH": "false"}
            run = make_run(manifest, root, cfg, "r")

            def fake_run(command, **kwargs):
                if command[0] == "shovill":
                    out = Path(command[command.index("--outdir") + 1]); out.mkdir(parents=True, exist_ok=True); (out / "contigs.fa").write_text(">c\n" + "A"*120 + "\n")
                elif Path(command[0]).name == "checkm2":
                    out = Path(command[command.index("--output-directory") + 1]); out.mkdir(parents=True, exist_ok=True)
                    kwargs["stdout"].write_text("predict started\n")
                    kwargs["stderr"].write_text("double free or corruption in CheckM2 metadata worker\n")
                    (out / "checkm2.log").write_text("internal CheckM2 worker crashed\n")
                    raise __import__("subprocess").CalledProcessError(134,command)
                elif command[0] == "prokka":
                    raise AssertionError("Prokka should not run after CheckM2 infrastructure failure")

            with patch("cleangene.workers.command_exists", return_value=True), patch("cleangene.workers.run", side_effect=fake_run):
                with self.assertRaisesRegex(SystemExit, "double free or corruption") as raised:
                    preprocess(run, 0)
            message=str(raised.exception)
            self.assertIn("predict started",message)
            self.assertIn("internal CheckM2 worker crashed",message)
            self.assertIn("command=",message)
            self.assertIn("exit_status=134",message)
            sample = run / "results" / "sample_data" / "iso1"
            self.assertFalse((sample / "qc.tsv").exists())
            self.assertFalse((run / "state" / "preprocess" / "iso1.done.json").exists())

    def test_checkm2_accepts_compressed_assembly_input_name(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); assembly=root/"contigs.fasta.gz"; db=root/EXPECTED_CHECKM2_DB_NAME; exe=_write_executable(root/"bin"/"checkm2")
            assembly.write_text("not really gzip but CheckM2 is mocked\n"); _write_db(db)
            from cleangene.workers import _run_checkm2
            commands=[]
            def fake_run(command, **kwargs):
                commands.append(command)
                out=Path(command[command.index("--output-directory")+1]); out.mkdir(parents=True,exist_ok=True)
                write_tsv(out/"quality_report.tsv",["Name","Completeness","Contamination"],[["iso1",99,1]])
            _run_root=root/"run"; logs=_run_root/"results"/"sample_data"/"iso1"/"logs"; logs.mkdir(parents=True)
            cfg={"CHECKM2_EXECUTABLE":str(exe),"CHECKM2_DB":str(db),"CHECKM2_VERSION":"CheckM2 version 1.1.0","CHECKM2_PREDICT_CPUS":"1","CHECKM2_MAX_INFLIGHT":"8"}
            with patch("cleangene.workers.run",side_effect=fake_run):
                self.assertEqual(_run_checkm2(assembly,_run_root/"results"/"sample_data"/"iso1"/"checkm2",logs,cfg,"iso1"),(99.0,1.0))
            self.assertTrue(( _run_root/"results"/"sample_data"/"iso1"/"checkm2"/"input"/"iso1.fasta.gz").is_symlink())

    def test_source_does_not_use_unsupported_hyphenated_cleanup_flag(self):
        root=Path(__file__).resolve().parents[1]
        offenders=[]
        for path in (root/"src").rglob("*.py"):
            if "--remove-intermediates" in path.read_text():
                offenders.append(str(path))
        self.assertEqual(offenders,[])


class RealCheckM2IntegrationTests(unittest.TestCase):
    def test_real_production_predict_smoke_when_companion_runtime_is_available(self):
        try:
            exe=resolve_checkm2_executable("")
            resolution=resolve_checkm2_db(DEFAULTS,allow_download=False)
            genome=bundled_test_genome(exe)
        except (ToolResolutionError, CheckM2DbError, CheckM2DbNotReady) as error:
            self.skipTest(f"real CheckM2 runtime/database unavailable: {error}")
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            input_path=checkm2_named_input_link(genome,root/"input","cleangene_checkm2_smoke")
            out=root/"predict"
            command=checkm2_predict_command(exe,input_path,out,resolution.path,1,checkm2_predict_capabilities(exe))
            env=__import__("os").environ.copy()
            for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","TF_NUM_INTRAOP_THREADS","TF_NUM_INTEROP_THREADS"):
                env[key]="1"
            completed=__import__("subprocess").run(command,capture_output=True,text=True,env=env,timeout=600)
            self.assertEqual(completed.returncode,0,completed.stderr[-2000:])
            completeness,contamination=parse_checkm2_quality_report(out/"quality_report.tsv","cleangene_checkm2_smoke")
            self.assertGreaterEqual(completeness,0)
            self.assertLessEqual(contamination,100)


if __name__ == "__main__":
    unittest.main()
