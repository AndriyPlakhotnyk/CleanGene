"""Resumable CheckM2 downloads, verified before publishing a usable database."""
from __future__ import annotations
import argparse
import hashlib
from http.client import HTTPException
import json
import re
import shutil
import subprocess
import tarfile
import time
from pathlib import Path
from urllib.request import Request, urlopen

from .checkm2 import EXPECTED_CHECKM2_DB_NAME, companion_python_for_checkm2
from .util import atomic_json, sha256


def compatible_database(executable: str) -> dict:
    script = '''
import json
from pathlib import Path
from checkm2.versionControl import VersionControl
from checkm2.defaultValues import DefaultValues
from checkm2.version import __version__
version, doi = VersionControl().return_highest_compatible_DB_version()
rows=json.loads((Path(DefaultValues.VERSION_PATH)/f'version_hashes_{__version__}.json').read_text())
row=next(r for r in rows if r['type']=='DIAMONDDB' and str(r['version'])==str(version))
print(json.dumps(row))
'''
    result = subprocess.run([str(companion_python_for_checkm2(executable)), '-c', script], check=True, capture_output=True, text=True, timeout=60)
    return json.loads(result.stdout)


def download_file(url: str, destination: Path, size: int, checksum: str, retries: int = 4) -> None:
    """Keep .part on connection failure and resume only a validated byte range."""
    algorithm, expected = checksum.split(':', 1)
    if algorithm not in {'md5','sha256'}: raise ValueError('Unsupported archive checksum')
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.part')
    def valid(path):
        digest = hashlib.new(algorithm)
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024*1024), b''): digest.update(block)
        return path.stat().st_size == size and digest.hexdigest() == expected
    if destination.is_file() and valid(destination): return
    for attempt in range(retries):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            if offset > size: partial.unlink(); offset = 0
            if offset < size:
                print(f'Download attempt {attempt+1}: resuming at {offset}/{size} bytes', flush=True)
                request = Request(url, headers={'Range': f'bytes={offset}-'} if offset else {})
                with urlopen(request, timeout=60) as response:
                    if offset and response.status == 206:
                        if not response.headers.get('Content-Range','').startswith(f'bytes {offset}-'):
                            raise ValueError('Server returned an incorrect resume range')
                        mode = 'ab'
                    elif response.status == 200: mode = 'wb'
                    else: raise ValueError(f'Unexpected download HTTP status: {response.status}')
                    with partial.open(mode) as handle:
                        shutil.copyfileobj(response, handle, length=1024*1024)
            if not valid(partial):
                if partial.stat().st_size >= size: partial.unlink()
                raise ValueError('Incomplete download or archive checksum mismatch')
            partial.replace(destination)
            return
        except (OSError, ValueError, HTTPException) as error:
            if attempt + 1 == retries: raise
            print(f'Download retry {attempt+1}/{retries}: {error}', flush=True)
            time.sleep(min(2**attempt, 8))


def download_database(executable: str, root: Path) -> Path:
    root = root.resolve(); root.mkdir(parents=True, exist_ok=True)
    specification = compatible_database(executable)
    match = re.fullmatch(r'10\.5281/zenodo\.(\d+)', specification['DOI'])
    if not match: raise ValueError('Installed CheckM2 does not specify a compatible Zenodo database')
    with urlopen(f'https://zenodo.org/api/records/{match[1]}', timeout=60) as response: metadata = json.load(response)
    archive_info = next(f for f in metadata['files'] if f['key']=='checkm2_database.tar.gz')
    cache = root / '.download'; cache.mkdir(exist_ok=True)
    archive = cache / 'checkm2_database.tar.gz'
    print(f"Downloading compatible database v{specification['version']} ({archive_info['size']} bytes)", flush=True)
    download_file(archive_info['links']['self'], archive, int(archive_info['size']), archive_info['checksum'])
    staging = cache / 'verified'; staging.mkdir(exist_ok=True)
    temporary = staging / EXPECTED_CHECKM2_DB_NAME
    print('Extracting and verifying database SHA-256', flush=True)
    with tarfile.open(archive, 'r:gz') as tar:
        members = [m for m in tar.getmembers() if Path(m.name).name == EXPECTED_CHECKM2_DB_NAME]
        if len(members)!=1 or not members[0].isfile(): raise ValueError('Archive lacks a unique regular database file')
        with tar.extractfile(members[0]) as source, temporary.open('wb') as target:
            shutil.copyfileobj(source, target, length=1024*1024)
    if sha256(temporary) != specification['sha256']:
        temporary.unlink()
        raise ValueError('Database checksum does not match installed CheckM2 version')
    destination = root / 'CheckM2_database' / EXPECTED_CHECKM2_DB_NAME
    destination.parent.mkdir(exist_ok=True)
    temporary.replace(destination)
    atomic_json(root / '.cleangene-download.json', {'status':'complete', 'database':str(destination), 'sha256':specification['sha256'], 'database_version':specification['version'], 'doi':specification['DOI'], 'archive_checksum':archive_info['checksum']})
    archive.unlink()
    print(f'Verified CheckM2 database ready: {destination}', flush=True)
    return destination


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',required=True); parser.add_argument('--path',required=True,type=Path)
    parser.add_argument('--no_write_json_db',action='store_true')
    args=parser.parse_args(); download_database(args.executable,args.path)

if __name__=='__main__': main()
