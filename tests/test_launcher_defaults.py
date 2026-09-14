import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from cleangene.cli import main
from cleangene.runtime import cleangene_project_root


class LauncherDefaultsTests(unittest.TestCase):
    def test_manifest_and_downsampling_only_use_checkout_and_local(self):
        with patch('cleangene.cli.run_command',return_value=0) as launch:
            self.assertEqual(main(['run','--manifest','input/example.tsv','--skip-downsampling']),0)
        args=launch.call_args.args[0]
        self.assertEqual(args.profile,'local')
        self.assertEqual(args.analysis_root,cleangene_project_root())
        self.assertTrue(args.skip_downsampling)

    def test_explicit_slurm_and_root_still_work(self):
        with patch('cleangene.cli.run_command',return_value=0) as launch:
            main(['run','--manifest','input/example.tsv','--profile','slurm','--analysis-root','/tmp/analysis'])
        args=launch.call_args.args[0]
        self.assertEqual(args.profile,'slurm')
        self.assertEqual(args.analysis_root,Path('/tmp/analysis'))

    def test_standalone_install_falls_back_to_current_directory(self):
        with patch('cleangene.cli.cleangene_project_root',return_value=None),patch('cleangene.cli.run_command',return_value=0) as launch:
            main(['run','--manifest','input/example.tsv'])
        self.assertEqual(launch.call_args.args[0].analysis_root,Path.cwd())

    def test_local_dry_run_does_not_launch_processing(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); manifest=root/'manifest.tsv'
            manifest.write_text('isolate_id\tR1\tR2\ni\t/r1.fastq.gz\t/r2.fastq.gz\n')
            with patch('cleangene.cli.local') as execute,patch('cleangene.cli.global_preflight') as preflight,contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(['run','--manifest',str(manifest),'--analysis-root',str(root),'--dry-run']),0)
            execute.assert_not_called();preflight.assert_not_called()
