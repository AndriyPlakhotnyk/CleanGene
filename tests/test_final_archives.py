import gzip
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_archive_discovery import archive_fixture
from cleangene.alignment_archive import alignment_digest, archive_own_alignment
from cleangene.final_archives import archive_bam, archive_pipeline_bams, restore_bam, restore_pipeline_bams
from cleangene.util import load_json, write_tsv


@unittest.skipUnless(shutil.which('samtools'), 'SAMtools required')
class FinalArchiveTests(unittest.TestCase):
    def test_all_archives_restore_and_leave_linked_inputs_and_utilities(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); ev,ref,bam=archive_fixture(root/'run/results/group')
            competitive=ev/'pangenome_reads.bam'; shutil.copy2(bam,competitive)
            expected=alignment_digest(competitive,None)
            external=root/'input.bam'; shutil.copy2(bam,external)
            (ev/'supplied.bam').symlink_to(external)
            utility=root/'run/results/utils/restored.bam'; utility.parent.mkdir(); shutil.copy2(bam,utility)
            inside=ev/"input.bam"; shutil.copy2(bam,inside)
            write_tsv(root/"run/provenance/manifest.tsv",["isolate_id","raw_bam"],[["i",str(inside)]])
            archive_own_alignment(ev)
            records=archive_pipeline_bams(root/'run')
            self.assertEqual(len(records),2)
            self.assertFalse(competitive.exists())
            self.assertTrue(external.is_file()); self.assertTrue(utility.is_file()); self.assertTrue(inside.is_file())
            self.assertEqual(len(archive_pipeline_bams(root/'run')),2)
            restored=restore_pipeline_bams(root/'run',root/'restored')
            self.assertEqual(len(restored),2)
            for row in restored:
                self.assertEqual(alignment_digest(Path(row['bam']),None),expected)
                self.assertTrue(Path(row['bam']+'.bai').is_file())

    def test_gzip_and_unsorted_roundtrip(self):
        for sort in ('coordinate','queryname'):
            with self.subTest(sort=sort),tempfile.TemporaryDirectory() as d:
                root=Path(d); ev,ref,bam=archive_fixture(root)
                if sort=='queryname':
                    sorted_bam=ev/'sorted.bam'
                    subprocess.run(['samtools','sort','-n','-o',str(sorted_bam),str(bam)],check=True)
                    sorted_bam.replace(bam)
                expected=alignment_digest(bam,None)
                zipped=Path(str(bam)+'.gz')
                with bam.open('rb') as source,gzip.open(zipped,'wb') as out: shutil.copyfileobj(source,out)
                bam.unlink()
                manifest=archive_bam(zipped)
                self.assertFalse(zipped.exists())
                self.assertEqual(bool(load_json(manifest)['index']),sort=='coordinate')
                restored=restore_bam(manifest,root/'restored.bam')
                self.assertEqual(alignment_digest(restored,None),expected)
                self.assertEqual(restore_bam(manifest,restored),restored)
                cram=manifest.parent/load_json(manifest)['cram']; cram.write_bytes(b'corrupt')
                with self.assertRaisesRegex(ValueError,'checksum mismatch'): restore_bam(manifest,root/'other.bam')

    def test_conversion_failure_preserves_bam(self):
        with tempfile.TemporaryDirectory() as d:
            ev,ref,bam=archive_fixture(Path(d)); before=bam.read_bytes()
            with patch('cleangene.final_archives.run',side_effect=subprocess.CalledProcessError(1,'samtools')):
                with self.assertRaises(subprocess.CalledProcessError): archive_bam(bam)
            self.assertEqual(bam.read_bytes(),before)
            self.assertTrue(Path(str(bam)+'.bai').is_file())


class DownsamplingCLITests(unittest.TestCase):
    def test_run_and_resume_flags_preserve_other_processing(self):
        from cleangene.cli import main, apply_cli_overrides
        for command,target in ((['run','--analysis-root','/tmp/example','--manifest','/tmp/example.tsv'],'run_command'),(['resume','--run-dir','/tmp/example'],'resume_command')):
            with self.subTest(command=command),patch('cleangene.cli.'+target,return_value=0) as execute:
                self.assertEqual(main(command+['--skip-downsampling']),0)
                cfg=apply_cli_overrides({'ASSEMBLER':'shovill','SKIP_TRIM':'false'},execute.call_args.args[0])
                self.assertEqual(cfg['SKIP_DOWNSAMPLING'],'true')
                self.assertEqual(cfg['ASSEMBLER'],'shovill')
                self.assertEqual(cfg['SKIP_TRIM'],'false')
