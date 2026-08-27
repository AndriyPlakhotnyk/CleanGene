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
from cleangene.runtime import assert_config_matches_runtime, cleangene_project_root
from cleangene.tools import ToolResolutionError, resolve_checkm2_executable, resolve_executable


def _write_executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/usr/bin/env sh\nprintf 'tool version 1.0\\n'\n")
    path.chmod(0o755)
    return path


class RuntimeResolutionTests(unittest.TestCase):
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
                with self.assertRaisesRegex(ToolResolutionError, "mamba env create --file environment.checkm2.yml"):
                    resolve_checkm2_executable(python_executable=python)

    def test_explicit_invalid_executable_fails_without_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            fallback = _write_executable(Path(d) / "bin" / "checkm2")
            with patch("cleangene.tools.shutil.which", return_value=str(fallback)):
                with self.assertRaisesRegex(ToolResolutionError, "Explicit checkm2 executable is invalid"):
                    resolve_executable("checkm2", str(Path(d) / "missing" / "checkm2"))

    def test_missing_checkm2_fails_before_controller_submission(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); manifest = root / "manifest.tsv"
            manifest.write_text("isolate_id\tgroup_id\tR1\tR2\ni1\tg\t/r1.fastq.gz\t/r2.fastq.gz\n")
            with patch("cleangene.cli.resolve_checkm2_executable", side_effect=ToolResolutionError("missing checkm2")), patch("cleangene.cli.slurm") as submit:
                with self.assertRaisesRegex(SystemExit, "before controller submission"):
                    with contextlib.redirect_stdout(StringIO()):
                        main(["run", "--manifest", str(manifest), "--analysis-root", str(root), "--run-id", "r", "--dry-run"])
            submit.assert_not_called()

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
