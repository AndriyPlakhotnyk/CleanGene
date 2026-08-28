import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from cleangene.runtime import cleangene_project_root


class InstallerTests(unittest.TestCase):
    def run_installer(
        self,
        environments: str,
        *arguments: str,
        existing_local: str | None = None,
        fail_doctor: bool = False,
    ):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "scripts").mkdir()
            (root / "config").mkdir()
            (root / "bin").mkdir()
            source = cleangene_project_root()
            shutil.copy(source / "scripts" / "install_or_update.sh", root / "scripts" / "install_or_update.sh")
            (root / "environment.yml").write_text("name: cleangene\n")
            (root / "environment.checkm2.yml").write_text("name: cleangene-checkm2\n")
            (root / "config" / "cleangene.arc.env").write_text("SLURM_ACCOUNT=\n")
            local = root / "config" / "cleangene.arc.local.env"
            if existing_local is not None:
                local.write_text(existing_local)
            log = root / "mamba.log"
            mamba = root / "bin" / "mamba"
            mamba.write_text(
                "#!/usr/bin/env bash\n"
                "printf '%s\\n' \"$*\" >> \"$MAMBA_LOG\"\n"
                "if [[ ${1:-} == env && ${2:-} == list ]]; then\n"
                "  for name in $MAMBA_ENVS; do printf '%s /tmp/%s\\n' \"$name\" \"$name\"; done\n"
                "fi\n"
                "if [[ ${MAMBA_FAIL_DOCTOR:-false} == true && $* == *'cleangene doctor'* ]]; then exit 17; fi\n"
            )
            mamba.chmod(0o755)
            env = {
                **os.environ,
                "PATH":f"{root / 'bin'}:{os.environ.get('PATH','')}",
                "MAMBA_LOG":str(log),
                "MAMBA_ENVS":environments,
                "MAMBA_FAIL_DOCTOR":"true" if fail_doctor else "false",
            }
            result = subprocess.run(["bash", "scripts/install_or_update.sh", *arguments], cwd=root, env=env, capture_output=True, text=True)
            return result, log.read_text().splitlines(), local.read_text()

    def test_create_mode_creates_both_environments_and_local_config(self):
        result, commands, local = self.run_installer("")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("env create -f environment.yml", commands)
        self.assertIn("env create -f environment.checkm2.yml", commands)
        self.assertFalse(any("env update" in command or "env remove" in command for command in commands))
        self.assertIn("run -n cleangene cleangene doctor --config config/cleangene.arc.local.env", commands)
        self.assertEqual(local, "SLURM_ACCOUNT=\n")

    def test_update_mode_updates_both_and_preserves_local_config(self):
        result, commands, local = self.run_installer("cleangene cleangene-checkm2", existing_local="PRIVATE=kept\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("env update -n cleangene -f environment.yml --prune", commands)
        self.assertIn("env update -n cleangene-checkm2 -f environment.checkm2.yml --prune", commands)
        self.assertFalse(any("env remove" in command for command in commands))
        self.assertEqual(local, "PRIVATE=kept\n")

    def test_recreate_mode_removes_and_creates_both(self):
        result, commands, _ = self.run_installer("cleangene cleangene-checkm2", "--recreate")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("env remove -n cleangene --yes", commands)
        self.assertIn("env remove -n cleangene-checkm2 --yes", commands)
        self.assertIn("env create -f environment.yml", commands)
        self.assertIn("env create -f environment.checkm2.yml", commands)
        self.assertFalse(any("env update" in command for command in commands))

    def test_doctor_failure_reports_installation_failure(self):
        result, commands, _ = self.run_installer(
            "cleangene cleangene-checkm2",
            fail_doctor=True,
        )
        self.assertEqual(result.returncode, 17)
        self.assertIn("run -n cleangene cleangene doctor --config config/cleangene.arc.local.env", commands)
        self.assertIn("ERROR: CleanGene installation/update failed.", result.stderr)


if __name__ == "__main__":
    unittest.main()
