"""Archive pipeline-generated BAMs without modifying linked or supplied inputs."""
from __future__ import annotations
import fcntl
import gzip
import os
import shutil
import subprocess
from pathlib import Path
from .alignment_archive import alignment_digest, sha256, restore_own_bam, verified_archive, archive_own_alignment
from .util import atomic_json, load_json, read_tsv, run, write_tsv


def archive_paths(bam: Path) -> tuple[Path, Path]:
    plain=bam.with_suffix('') if bam.name.lower().endswith('.bam.gz') else bam
    cram=plain.with_suffix('.cram')
    return cram, cram.with_suffix('.cram.archive.json')


def archive_bam(bam: Path, threads: int = 1) -> Path:
    """Store sequence bases in CRAM when the alignment reference is unknown.

    This also supports name-collated intermediates: their order is preserved and
    an index is created only for coordinate-sorted alignments.
    """
    cram, manifest=archive_paths(bam)
    with manifest.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        source=bam
        unpacked=bam.with_name('.'+bam.name+'.archive-input.bam')
        try:
            if bam.name.lower().endswith('.bam.gz'):
                with gzip.open(bam,'rb') as src, unpacked.open('wb') as dst: shutil.copyfileobj(src,dst)
                source=unpacked
            expected,count=alignment_digest(source,None)
            if manifest.is_file():
                previous=load_json(manifest)
                if previous['source_name'] != bam.name and (previous['alignment_sha256'],previous['alignment_count']) != (expected,count):
                    raise ValueError(f'Conflicting BAM archive names: {bam} and {previous["source_name"]}')
            header=subprocess.check_output(['samtools','view','-H',str(source)],text=True)
            indexed=any(line.startswith('@HD\t') and 'SO:coordinate' in line.split('\t') for line in header.splitlines())
            temporary=cram.with_name('.'+cram.name+'.tmp.cram')
            run(['samtools','view','--no-PG','-C','--output-fmt-option','no_ref=1',
                 '--output-fmt-option','store_md=1','--output-fmt-option','store_nm=1',
                 '-@',str(threads),'-o',str(temporary),str(source)])
            run(['samtools','quickcheck','-u',str(temporary)])
            if alignment_digest(temporary,None)!=(expected,count): raise ValueError(f'CRAM verification failed; retaining {bam}')
            if indexed: run(['samtools','index',str(temporary)])
            temporary.replace(cram)
            index=Path(str(cram)+'.crai')
            if indexed: Path(str(temporary)+'.crai').replace(index)
            else: index.unlink(missing_ok=True)
            atomic_json(manifest,{'version':1,'cram':cram.name,'cram_sha256':sha256(cram),
                        'index':index.name if indexed else '', 'index_sha256':sha256(index) if indexed else '',
                        'reference_mode':'stored_bases','source_name':bam.name,'bam_name':cram.with_suffix('.bam').name,
                        'alignment_sha256':expected,'alignment_count':count,'lossless_verified':True})
            bam.unlink()
            for suffix in ('.bai','.csi'):
                Path(str(bam)+suffix).unlink(missing_ok=True)
                if bam.name.lower().endswith('.bam.gz'): Path(str(bam.with_suffix(''))+suffix).unlink(missing_ok=True)
            return manifest
        finally: unpacked.unlink(missing_ok=True)


def restore_bam(manifest: Path, output: Path, threads: int = 1) -> Path:
    """Restore and verify any BAM in a final-archive manifest."""
    meta=load_json(manifest); cram=manifest.parent/meta['cram']
    if sha256(cram)!=meta['cram_sha256']: raise ValueError(f'CRAM checksum mismatch: {cram}')
    if meta['index'] and sha256(manifest.parent/meta['index'])!=meta['index_sha256']: raise ValueError('CRAI checksum mismatch')
    expected=(meta['alignment_sha256'],meta['alignment_count'])
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.with_suffix('.restore.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if output.exists():
            if alignment_digest(output,None)!=expected: raise ValueError(f'Refusing to overwrite different BAM: {output}')
        else:
            temporary=output.with_name('.'+output.name+'.tmp.bam')
            run(['samtools','view','--no-PG','-b','--input-fmt-option','decode_md=0','-@',str(threads),'-o',str(temporary),str(cram)])
            if alignment_digest(temporary,None)!=expected: raise ValueError('Restored BAM failed verification')
            temporary.replace(output)
        if meta['index']: run(['samtools','index',str(output)])
    return output


def pipeline_files(run_dir: Path):
    root=(run_dir/'results').resolve()
    manifest=run_dir/'provenance/manifest.tsv'
    supplied={Path(row['raw_bam']).expanduser().resolve() for row in read_tsv(manifest) if row.get('raw_bam')} if manifest.is_file() else set()
    for directory,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=[d for d in dirs if not Path(directory,d).is_symlink() and not (Path(directory)==root and d=='utils')]
        for name in files:
            path=Path(directory,name)
            if not path.is_symlink() and path.resolve() not in supplied: yield path


def archive_pipeline_bams(run_dir: Path, threads: int = 1) -> list[dict]:
    for path in list(pipeline_files(run_dir)):
        if path.name=='own_assembly_reads.bam' and (path.parent/'own_assembly_reference.fasta').is_file():
            archive_own_alignment(path.parent,threads)
        elif path.name.lower().endswith(('.bam','.bam.gz')): archive_bam(path,threads)
    records=[]
    for path in pipeline_files(run_dir):
        if path.name.endswith('.cram.archive.json'):
            meta=load_json(path)
            records.append({'manifest':str(path.relative_to(run_dir.resolve())), 'cram':str((path.parent/meta['cram']).relative_to(run_dir.resolve())), 'reference_mode':meta['reference_mode'],'alignment_count':meta['alignment_count']})
        elif path.name=='own_assembly_reads.archive.json':
            cram,reference,meta=verified_archive(path.parent)
            records.append({'manifest':str(path.relative_to(run_dir.resolve())),'cram':str(cram.relative_to(run_dir.resolve())),'reference_mode':'retained_assembly','alignment_count':meta['alignment_count']})
    records.sort(key=lambda row:row['manifest'])
    write_tsv(run_dir/'results/cohort/alignment_archives.tsv',['manifest','cram','reference_mode','alignment_count'],records)
    return records


def restore_pipeline_bams(run_dir: Path, output: Path, threads: int = 1) -> list[dict]:
    records=[]
    for path in pipeline_files(run_dir):
        relative=path.parent.relative_to((run_dir/'results').resolve())
        if path.name.endswith('.cram.archive.json'):
            meta=load_json(path); bam=restore_bam(path,output/relative/meta['bam_name'],threads)
        elif path.name=='own_assembly_reads.archive.json':
            bam=restore_own_bam(path.parent,output/relative/'own_assembly_reads.bam',threads)
        else: continue
        records.append({'archive_manifest':str(path),'bam':str(bam)})
    write_tsv(output/'restored_bams.tsv',['archive_manifest','bam'],records)
    return records
