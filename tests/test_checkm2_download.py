import hashlib
import io
import json
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from cleangene.checkm2_download import download_file, download_database
from cleangene.checkm2 import EXPECTED_CHECKM2_DB_NAME, find_managed_checkm2_db


class Response(io.BytesIO):
    def __init__(self, body, status=200, headers=None):
        super().__init__(body); self.status=status; self.headers=headers or {}


class DownloadTests(unittest.TestCase):
    def test_resume_and_checksum_verified_reuse(self):
        body=b'correct archive contents'; checksum='md5:'+hashlib.md5(body).hexdigest()
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'archive'; path.with_name('archive.part').write_bytes(body[:7])
            with patch('cleangene.checkm2_download.urlopen',return_value=Response(body[7:],206,{'Content-Range':f'bytes 7-{len(body)-1}/{len(body)}'})) as fetch:
                download_file('https://example.test/data',path,len(body),checksum)
                self.assertEqual(fetch.call_args.args[0].get_header('Range'),'bytes=7-')
                self.assertEqual(path.read_bytes(),body)
            with patch('cleangene.checkm2_download.urlopen',side_effect=AssertionError('unnecessary download')):
                download_file('https://example.test/data',path,len(body),checksum)

    def test_server_ignoring_range_restarts_instead_of_appending(self):
        body=b'abcdefgh'
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'archive'; path.with_name('archive.part').write_bytes(body[:2])
            with patch('cleangene.checkm2_download.urlopen',return_value=Response(body)):
                download_file('https://example.test/data',path,len(body),'sha256:'+hashlib.sha256(body).hexdigest())
            self.assertEqual(path.read_bytes(),body)

    def test_bad_checksum_never_publishes_and_incomplete_bytes_remain_resumable(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'archive'
            with patch('cleangene.checkm2_download.urlopen',return_value=Response(b'bad')):
                with self.assertRaisesRegex(ValueError,'checksum'):
                    download_file('https://example.test/data',path,5,'md5:'+hashlib.md5(b'right').hexdigest(),retries=1)
            self.assertFalse(path.exists())
            self.assertEqual(path.with_name('archive.part').read_bytes(),b'bad')

    def test_database_publication_checks_installed_version_hash(self):
        body=b'verified diamond database'
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); archive=root/'source.tar.gz'
            with tarfile.open(archive,'w:gz') as tar:
                info=tarfile.TarInfo('CheckM2_database/'+EXPECTED_CHECKM2_DB_NAME); info.size=len(body)
                tar.addfile(info,io.BytesIO(body))
            specification={'DOI':'10.5281/zenodo.14897628','version':'3','sha256':hashlib.sha256(body).hexdigest()}
            metadata={'files':[{'key':'checkm2_database.tar.gz','size':archive.stat().st_size,'checksum':'md5:test','links':{'self':'https://example.test/archive'}}]}
            def copy(url,destination,*args): shutil.copyfile(archive,destination)
            def response(*args,**kwargs): return Response(json.dumps(metadata).encode())
            with patch('cleangene.checkm2_download.compatible_database',return_value=specification), patch('cleangene.checkm2_download.urlopen',side_effect=response), patch('cleangene.checkm2_download.download_file',side_effect=copy):
                db=download_database('checkm2',root/'managed')
                self.assertEqual(db.read_bytes(),body)
                specification['sha256']='wrong'
                with self.assertRaisesRegex(ValueError,'installed CheckM2'): download_database('checkm2',root/'managed')
                self.assertEqual(db.read_bytes(),body)

    def test_unverified_staging_is_never_a_managed_database(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); stage=root/'.download/verified'; stage.mkdir(parents=True)
            (stage/EXPECTED_CHECKM2_DB_NAME).write_bytes(b'partial extraction')
            self.assertIsNone(find_managed_checkm2_db(root))

class CompanionEnvironmentTests(unittest.TestCase):
    def test_companion_tools_precede_main_environment_without_mutating_parent(self):
        import os
        from cleangene.checkm2 import checkm2_subprocess_environment
        with tempfile.TemporaryDirectory() as d:
            exe=Path(d)/'companion/bin/checkm2'
            with patch.dict(os.environ,{'PATH':'/main/bin','OMP_NUM_THREADS':'64'}):
                env=checkm2_subprocess_environment(exe)
                self.assertEqual(env['PATH'],str(exe.parent)+os.pathsep+'/main/bin')
                self.assertEqual(env['OMP_NUM_THREADS'],'1')
                self.assertEqual(os.environ['PATH'],'/main/bin')
                self.assertEqual(os.environ['OMP_NUM_THREADS'],'64')

    def test_bundled_tst_genome_is_found_by_real_python_import(self):
        import os
        import sys
        from cleangene.checkm2 import bundled_test_genome
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); package=root/'modules/checkm2'; (package/'testrun').mkdir(parents=True)
            (package/'__init__.py').touch(); genome=package/'testrun/TEST2.tst'; genome.write_text('>genome\nACGT\n')
            exe=root/'bin/checkm2'; exe.parent.mkdir(); exe.touch(); (exe.parent/'python').symlink_to(sys.executable)
            with patch.dict(os.environ,{'PYTHONPATH':str(package.parent)}):
                self.assertEqual(bundled_test_genome(exe),genome)
