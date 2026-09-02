import contextlib
import os
import subprocess
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from cleangene.checkm2 import checkm2_database_root
from cleangene.cli import main
from cleangene.kraken import managed_kraken2_db_path
from cleangene.runtime import assert_config_matches_runtime, cleangene_project_root, record_runtime_provenance, verify_worker_runtime
from cleangene.tools import (
    CHECKM2_VERSION_TIMEOUT_SECONDS,
    ToolResolutionError,
    executable_version,
    resolve_checkm2_executable,
    resolve_executable,
)


def _write_executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/usr/bin/env sh\nif [ \"$1 $2\" = 'predict --help' ]; then printf '%s\\n' '--input --output-directory --database_path --threads --remove_intermediates --force'; else printf 'tool version 1.0\\n'; fi\n")
    path.chmod(0o755)
    return path


class RuntimeResolutionTests(unittest.TestCase):
    def test_worker_runtime_mismatch_is_detected(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d)/"run"; (run/"provenance").mkdir(parents=True)
            record_runtime_provenance(run,{})
            data=__import__("cleangene.util",fromlist=["load_json"]).load_json(run/"provenance"/"runtime.json")
            data["workers_py_sha256"]="old-worker"
            __import__("cleangene.util",fromlist=["atomic_json"]).atomic_json(run/"provenance"/"runtime.json",data)
            with self.assertRaisesRegex(SystemExit,"CleanGene runtime mismatch"):
                verify_worker_runtime(run)

    def test_checkm2_version_probe_allows_slow_hpc_startup(self):
        completed = subprocess.CompletedProcess(["checkm2", "--version"], 0, "1.1.0\n", "")
        with patch("cleangene.tools.subprocess.run", return_value=completed) as run:
            self.assertEqual(executable_version("checkm2", "CheckM2"), "1.1.0")
        run.assert_called_once_with(
            ["checkm2", "--version"],
            capture_output=True,
            text=True,
            timeout=CHECKM2_VERSION_TIMEOUT_SECONDS,
        )

    def test_checkm2_version_probe_reports_timeout(self):
        with patch(
            "cleangene.tools.subprocess.run",
            side_effect=subprocess.TimeoutExpired(["checkm2", "--version"], CHECKM2_VERSION_TIMEOUT_SECONDS),
        ):
            with self.assertRaisesRegex(ToolResolutionError, "within 300 seconds"):
                executable_version("checkm2", "CheckM2")

    def test_resolve_executable_from_path(self):
        with tempfile.TemporaryDirectory() as d:
            exe = _write_executable(Path(d) / "bin" / "checkm2")
            with patch("cleangene.tools.shutil.which", return_value=str(exe)):
                self.assertEqual(resolve_executable("checkm2"), exe.resolve())

    def test_resolve_executable_beside_python(self):
        with tempfile.TemporaryDirectory() as d:
            bindir = Path(d) / "env" / "bin"; python = _write_executable(bindir / "python"); exe = _write_executable(bindir / "checkm2")
            with patch("cleangene.tools.shutil.which", return_value=None):
                self.assertEqual(resolve_executable("checkm2", python_executable=python), exe.resolve())

    def test_resolve_checkm2_from_companion_environment(self):
        with tempfile.TemporaryDirectory() as d:
            envs = Path(d) / "envs"
            python = _write_executable(envs / "cleangene" / "bin" / "python")
            exe = _write_executable(envs / "cleangene-checkm2" / "bin" / "checkm2")
            with patch("cleangene.tools.shutil.which", return_value=None):
                self.assertEqual(resolve_checkm2_executable(python_executable=python), exe.resolve())

    def test_missing_checkm2_names_companion_environment_repair(self):
        with tempfile.TemporaryDirectory() as d:
            python = _write_executable(Path(d) / "envs" / "cleangene" / "bin" / "python")
            with patch("cleangene.tools.shutil.which", return_value=None):
                with self.assertRaisesRegex(ToolResolutionError, "bash scripts/install_or_update.sh --recreate"):
                    resolve_checkm2_executable(python_executable=python)

    def test_explicit_checkm2_override_wins(self):
        with tempfile.TemporaryDirectory() as d:
            explicit = _write_executable(Path(d) / "custom" / "checkm2")
            fallback = _write_executable(Path(d) / "bin" / "checkm2")
            with patch("cleangene.tools.shutil.which", return_value=str(fallback)):
                self.assertEqual(resolve_checkm2_executable(str(explicit)), explicit.resolve())

    def test_explicit_invalid_executable_fails_without_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            fallback = _write_executable(Path(d) / "bin" / "checkm2")
            with patch("cleangene.tools.shutil.which", return_value=str(fallback)):
                with self.assertRaisesRegex(ToolResolutionError, "Explicit checkm2 executable is invalid"):
                    resolve_executable("checkm2", str(Path(d) / "missing" / "checkm2"))

    def test_launcher_does_not_resolve_checkm2_before_controller_submission(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); manifest = root / "manifest.tsv"
            manifest.write_text("isolate_id\tgroup_id\tR1\tR2\ni1\tg\t/r1.fastq.gz\t/r2.fastq.gz\n")
            with patch("cleangene.cli.resolve_checkm2_executable", side_effect=ToolResolutionError("missing checkm2")), patch("cleangene.cli.slurm") as submit:
                with contextlib.redirect_stdout(StringIO()):
                    self.assertEqual(main(["run", "--manifest", str(manifest), "--analysis-root", str(root), "--run-id", "r", "--dry-run"]),0)
            submit.assert_called_once()

    def test_doctor_reports_automatic_missing_checkm2_database(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            exe = _write_executable(root / "bin" / "checkm2")
            config = root / "arc.local.env"
            config.write_text(
                f"CHECKM2_EXECUTABLE={exe}\n"
                f"CHECKM2_DATABASE_ROOT={root / 'checkm2'}\n"
                "CHECKM2_AUTO_DOWNLOAD=true\n"
                "TAXONOMY_MODE=off\n"
            )
            output = StringIO()
            environment = {**os.environ, "CONDA_DEFAULT_ENV":"cleangene", "CONDA_PREFIX":str(root / "envs" / "cleangene")}
            with patch.dict(os.environ, environment, clear=True), patch("cleangene.cli.sys.executable", str(root / "envs" / "cleangene" / "bin" / "python")), patch("cleangene.cli.command_exists", return_value=True), contextlib.redirect_stdout(output):
                result = main(["doctor", "--config", str(config)])
            self.assertEqual(result, 0, output.getvalue())
            self.assertIn("CheckM2 executable: READY", output.getvalue())
            self.assertIn("CheckM2 database: not present; will be created by checkm2_db_setup", output.getvalue())

    def test_doctor_checkm2_failure_names_repair_command(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            config = root / "arc.local.env"
            config.write_text("CHECKM2_MODE=required\nTAXONOMY_MODE=off\n")
            output = StringIO()
            environment = {**os.environ, "CONDA_DEFAULT_ENV":"cleangene", "CONDA_PREFIX":str(root / "envs" / "cleangene")}
            with patch.dict(os.environ, environment, clear=True), patch("cleangene.cli.sys.executable", str(root / "envs" / "cleangene" / "bin" / "python")), patch("cleangene.cli.command_exists", return_value=True), patch("cleangene.cli.resolve_checkm2_executable", side_effect=ToolResolutionError("missing checkm2")), contextlib.redirect_stdout(output):
                result = main(["doctor", "--config", str(config)])
            self.assertEqual(result, 2)
            self.assertIn("Fix: bash scripts/install_or_update.sh --recreate", output.getvalue())

    def test_doctor_deep_checkm2_refuses_login_node_when_slurm_available(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            exe=_write_executable(root/"bin"/"checkm2")
            db=root/"checkm2"/"CheckM2_database"/"uniref100.KO.1.dmnd"
            db.parent.mkdir(parents=True); db.write_text("db\n")
            config=root/"arc.local.env"
            config.write_text(
                f"CHECKM2_MODE=required\nCHECKM2_EXECUTABLE={exe}\nCHECKM2_DATABASE_ROOT={root/'checkm2'}\n"
                "CHECKM2_CPUS=16\nCHECKM2_MEM=128G\nCHECKM2_TIME=24:00:00\nTAXONOMY_MODE=off\n"
            )
            output=StringIO()
            environment={**os.environ,"CONDA_DEFAULT_ENV":"cleangene","CONDA_PREFIX":str(root/"envs"/"cleangene")}
            environment.pop("SLURM_JOB_ID",None)
            with patch.dict(os.environ,environment,clear=True), \
                 patch("cleangene.cli.sys.executable",str(root/"envs"/"cleangene"/"bin"/"python")), \
                 patch("cleangene.cli.command_exists",return_value=True), \
                 contextlib.redirect_stdout(output):
                result=main(["doctor","--config",str(config),"--deep-checkm2"])
            text=output.getvalue()
            self.assertEqual(result,2)
            self.assertIn("CheckM2 database: READY",text)
            self.assertIn("CheckM2 deep runtime verification: ERROR",text)
            self.assertIn("salloc --cpus-per-task=16 --mem=128G --time=24:00:00",text)

    def test_config_from_other_checkout_fails_without_external_database_root(self):
        active = cleangene_project_root()
        self.assertIsNotNone(active)
        with tempfile.TemporaryDirectory() as d:
            other = Path(d) / "CleanGeneOther"; (other / "src" / "cleangene").mkdir(parents=True); (other / "pyproject.toml").write_text("[project]\nname='fake'\n")
            config = other / "config" / "cleangene.arc.env"; config.parent.mkdir(); config.write_text("CHECKM2_MODE=required\n")
            with self.assertRaisesRegex(SystemExit, "checkout mismatch"):
                assert_config_matches_runtime(config, {"CHECKM2_MODE": "required"})
            external = Path(d) / "shared-databases"
            assert_config_matches_runtime(config, {"CLEANGENE_DATABASE_ROOT": str(external)})

    def test_generic_database_root_feeds_kraken_and_checkm2(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "managed"
            kraken, size, _ = managed_kraken2_db_path({"CLEANGENE_DATABASE_ROOT": str(root), "KRAKEN2_DATABASE_SIZE": "standard-8"})
            self.assertEqual(size, "standard-8")
            self.assertEqual(kraken, (root / "kraken2_standard-8").resolve())
            self.assertEqual(checkm2_database_root({"CLEANGENE_DATABASE_ROOT": str(root)}), (root / "checkm2").resolve())


class RepositoryPrivacyTests(unittest.TestCase):
    def test_private_runtime_paths_are_gitignored(self):
        root = cleangene_project_root()
        for path in ("config/test.local.env", "config/site.private.env", "databases/example", "runs/example"):
            result = subprocess.run(["git", "check-ignore", path], cwd=root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, path)

    def test_tracked_public_config_has_no_private_site_values(self):
        root = cleangene_project_root()
        config = root / "config" / "cleangene.arc.env"
        text = config.read_text()
        self.assertIn('SLURM_ACCOUNT=""', text)
        self.assertIn('CHECKM2_EXECUTABLE=""', text)
        self.assertIn('CHECKM2_DB=""', text)
        self.assertIn('CLEANGENE_DATABASE_ROOT=""', text)
        forbidden = tuple(f"/{name}/" for name in ("home", "work", "scratch"))
        for value in forbidden:
            self.assertNotIn(value, text)

    def test_arc_config_has_unique_keys_and_unchanged_qc_thresholds(self):
        root = cleangene_project_root()
        lines = (root / "config" / "cleangene.arc.env").read_text().splitlines()
        keys = [line.split("=", 1)[0].strip() for line in lines if line.strip() and not line.lstrip().startswith("#") and "=" in line]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertNotIn("CHECKM2_BATCH_SIZE", keys)
        self.assertIn("CHECKM2_PREDICT_CPUS", keys)
        self.assertIn("CHECKM2_MAX_INFLIGHT", keys)
        expected = {
            "QC_MAX_CONTIGS_PASS":"300", "QC_MAX_CONTIGS_FAIL":"1000",
            "QC_MIN_N50_PASS":"25000", "QC_MIN_N50_FAIL":"5000",
            "QC_MIN_COVERAGE_PASS":"20", "QC_MIN_COVERAGE_FAIL":"10",
            "QC_MIN_READ_LENGTH_PASS":"120", "QC_MIN_READ_LENGTH_FAIL":"",
            "QC_MIN_MEAN_BASE_QUALITY_PASS":"30", "QC_MIN_MEAN_BASE_QUALITY_FAIL":"",
            "QC_MIN_COMPLETENESS_PASS":"90", "QC_MIN_COMPLETENESS_FAIL":"80",
            "QC_MAX_CHECKM2_CONTAMINATION_PASS":"5", "QC_MAX_CHECKM2_CONTAMINATION_FAIL":"10",
            "QC_MAX_KRAKEN_CONTAMINATION_FAIL":"5", "QC_PROFILE_FILE":"",
        }
        from cleangene.config import read_env
        from cleangene.cli import _doctor_config_errors
        cfg = read_env(root / "config" / "cleangene.arc.env")
        self.assertEqual({key:cfg[key] for key in expected}, expected)
        self.assertEqual(_doctor_config_errors(cfg), [])

    def test_tracked_text_files_do_not_contain_obvious_private_paths(self):
        root = cleangene_project_root()
        result = subprocess.run(["git", "ls-files"], cwd=root, capture_output=True, text=True, check=True)
        suffixes = {".py", ".md", ".txt", ".yml", ".yaml", ".env", ".sh", ".toml", ".gitignore"}
        offenders = []
        for rel in result.stdout.splitlines():
            if rel == "tests/test_runtime.py":
                continue
            path = root / rel
            if path.suffix not in suffixes and path.name != ".gitignore":
                continue
            text = path.read_text(errors="ignore")
            for marker in tuple(f"/{name}/" for name in ("home", "work", "scratch")):
                if marker in text:
                    offenders.append(f"{rel}:{marker}")
        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
