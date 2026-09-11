"""Lossless, reference-backed archival of sample-assembly alignments."""
from __future__ import annotations

import fcntl
import gzip
import hashlib
import shutil
import subprocess
from pathlib import Path

from .util import atomic_json, load_json, run

ARCHIVE_VERSION = 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''): digest.update(block)
    return digest.hexdigest()


def retain_reference(source: Path, evidence: Path) -> Path:
    """Keep the exact uncompressed FASTA, independently of assembly cleanup."""
    evidence.mkdir(parents=True, exist_ok=True)
    target = evidence / 'own_assembly_reference.fasta'
    if source.resolve() == target.resolve(): return target
    temporary = target.with_suffix('.fasta.tmp')
    opener = gzip.open if source.suffix == '.gz' else open
    with opener(source, 'rb') as src, temporary.open('wb') as dst: shutil.copyfileobj(src, dst)
    if target.exists() and sha256(target) == sha256(temporary): temporary.unlink()
    else:
        if (evidence / "own_assembly_reads.archive.json").exists():
            temporary.unlink()
            raise ValueError("Assembly differs from retained CRAM reference; use a new run/evidence directory")
        temporary.replace(target)
        for suffix in ('.fai', '.amb', '.ann', '.bwt', '.pac', '.sa'): Path(str(target) + suffix).unlink(missing_ok=True)
    atomic_json(evidence / 'own_assembly_reference.json', {'reference': target.name, 'sha256': sha256(target), 'source_sha256': sha256(source), 'source_name': source.name})
    return target


def alignment_digest(path: Path, reference: Path) -> tuple[str, int]:
    """Compare every alignment, quality and optional tag, ignoring tag ordering."""
    command = ['samtools', 'view', '--no-PG', '-T', str(reference)]
    if path.suffix == '.cram': command += ['--input-fmt-option', 'decode_md=0']
    command += [str(path)]
    digest = hashlib.sha256(); count = 0
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.rstrip('\n').split('\t')
            digest.update(('\t'.join([*fields[:11], *sorted(fields[11:])]) + '\n').encode()); count += 1
        status = process.wait()
    if status: raise subprocess.CalledProcessError(status, command)
    return digest.hexdigest(), count


def verified_archive(evidence: Path) -> tuple[Path, Path, dict]:
    metadata = load_json(evidence / 'own_assembly_reads.archive.json')
    cram = evidence / metadata['cram']; reference = evidence / metadata['reference']; index = evidence / metadata['index']
    for path, key in ((reference, 'reference_sha256'), (cram, 'cram_sha256'), (index, 'index_sha256')):
        if not path.is_file() or sha256(path) != metadata[key]: raise ValueError(f'Archived evidence checksum mismatch: {path}')
    return cram, reference, metadata


def archive_own_alignment(evidence: Path, threads: int = 1) -> dict:
    """Publish verified CRAM/CRAI/FASTA before removing the replaceable BAM/BAI."""
    evidence.mkdir(parents=True, exist_ok=True)
    with (evidence / '.archive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        bam = evidence / 'own_assembly_reads.bam'
        if not bam.is_file(): return verified_archive(evidence)[2]
        reference = evidence / 'own_assembly_reference.fasta'
        if not reference.is_file(): raise FileNotFoundError(f'Exact alignment reference missing: {reference}')
        run(['samtools', 'faidx', str(reference)])
        expected, count = alignment_digest(bam, reference)
        cram = evidence / 'own_assembly_reads.cram'; temporary = evidence / '.own_assembly_reads.tmp.cram'
        index = Path(str(cram) + '.crai'); tmp_index = Path(str(temporary) + '.crai')
        run(['samtools', 'view', '--no-PG', '-C', '-T', str(reference), '-@', str(max(1, threads)),
             '--output-fmt-option', 'store_md=1', '--output-fmt-option', 'store_nm=1', '-o', str(temporary), str(bam)])
        run(['samtools', 'index', str(temporary)])
        run(['samtools', 'quickcheck', '-v', str(temporary)])
        observed, observed_count = alignment_digest(temporary, reference)
        if (observed, observed_count) != (expected, count): raise ValueError('CRAM round-trip verification failed; original BAM retained')
        temporary.replace(cram); tmp_index.replace(index)
        metadata = {'version': ARCHIVE_VERSION, 'cram': cram.name, 'index': index.name, 'reference': reference.name,
                    'reference_sha256': sha256(reference), 'cram_sha256': sha256(cram), 'index_sha256': sha256(index),
                    'alignment_sha256': expected, 'alignment_count': count, 'lossless_verified': True}
        atomic_json(evidence / 'own_assembly_reads.archive.json', metadata)
        bam.unlink()
        Path(str(bam) + '.bai').unlink(missing_ok=True)
        return metadata


def restore_own_bam(evidence: Path, output: Path, threads: int = 1) -> Path:
    cram, reference, metadata = verified_archive(evidence)
    output.parent.mkdir(parents=True, exist_ok=True)
    with (output.parent / ('.' + output.name + '.restore.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if output.is_file():
            if alignment_digest(output, reference) != (metadata['alignment_sha256'], metadata['alignment_count']):
                raise ValueError(f'Refusing to overwrite different BAM: {output}')
            if not Path(str(output) + '.bai').is_file(): run(['samtools', 'index', str(output)])
            return output
        temporary = output.with_name('.' + output.name + '.tmp.bam')
        run(['samtools', 'view', '--no-PG', '-b', '-T', str(reference), '--input-fmt-option', 'decode_md=0',
             '-@', str(max(1, threads)), '-o', str(temporary), str(cram)])
        run(['samtools', 'index', str(temporary)])
        if alignment_digest(temporary, reference) != (metadata['alignment_sha256'], metadata['alignment_count']):
            raise ValueError('Restored BAM failed alignment verification')
        temporary.replace(output); Path(str(temporary) + '.bai').replace(Path(str(output) + '.bai'))
    return output
