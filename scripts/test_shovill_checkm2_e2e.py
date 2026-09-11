#!/usr/bin/env python3
"""Real two-isolate Shovill/CheckM2 pipeline and resume test (requires installed tools)."""
from __future__ import annotations
import argparse
import gzip
import json
import random
import subprocess
import sys
from pathlib import Path

from cleangene.alignment_archive import verified_archive
from cleangene.checkm2 import bundled_test_genome
from cleangene.fasta import read_fasta
from cleangene.util import read_tsv, sha256, write_tsv


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',required=True,type=Path)
    parser.add_argument('--checkm2-executable',required=True,type=Path)
    parser.add_argument('--database-root',required=True,type=Path)
    parser.add_argument('--genome',type=Path)
    parser.add_argument('--threads',type=int,default=2)
    args=parser.parse_args()
    root=args.work_dir.resolve(); root.mkdir(parents=True,exist_ok=True)
    if (root/'runs/e2e').exists(): raise SystemExit('Use a fresh work directory; existing evidence is preserved')
    genome=args.genome or bundled_test_genome(args.checkm2_executable)
    sequence=max(read_fasta(genome).values(),key=len)[:200000]
    if len(sequence)<50000: raise SystemExit('Test reference must have a contig of at least 50 kb')
    manifest=[]
    for isolate,seed in (('sample_a',81),('sample_b',82)):
        rng=random.Random(seed); paths=[root/f'{isolate}_R{mate}.fastq.gz' for mate in (1,2)]
        with gzip.open(paths[0],'wt') as a,gzip.open(paths[1],'wt') as b:
            for i in range(len(sequence)*200//300):
                start=rng.randrange(len(sequence)-400)
                pair=(sequence[start:start+150],sequence[start+250:start+400].translate(str.maketrans('ACGT','TGCA'))[::-1])
                for handle,read in zip((a,b),pair): handle.write(f'@r{i}\n{read}\n+\n'+('I'*len(read))+'\n')
        manifest.append([isolate,'e2e',*map(str,paths)])
    manifest_path=root/'manifest.tsv'; write_tsv(manifest_path,['isolate_id','group_id','R1','R2'],manifest)
    # This cropped-genome fixture tests software execution, not completeness accuracy.
    # Relax only the completeness exclusion gate, retaining measured CheckM2 values.
    cfg={'TAXONOMY_MODE':'off','CHECKM2_MODE':'required','CHECKM2_EXECUTABLE':str(args.checkm2_executable.resolve()),
         'CHECKM2_DATABASE_ROOT':str(args.database_root.resolve()),'CPUS':str(args.threads),
         'CHECKM2_PREDICT_CPUS':str(args.threads),'SHOVILL_MEMORY_GB':'8',
         'VALIDATION_CPUS':str(args.threads),'ARBITRATION_CPUS':str(args.threads),
         'PANAROO_SMALL_CPUS':str(args.threads),'PREPROCESS_USE_NODE_LOCAL_SCRATCH':'false',
         'QC_MIN_COMPLETENESS_PASS':'0','QC_MIN_COMPLETENESS_FAIL':'0'}
    config=root/'config.env'; config.write_text(''.join(f'{key}="{value}"\n' for key,value in cfg.items()))
    command=[sys.executable,'-m','cleangene','run','--profile','local','--manifest',str(manifest_path),
             '--config',str(config),'--analysis-root',str(root),'--run-id','e2e','--assembler','shovill']
    subprocess.run(command,check=True)
    run=root/'runs/e2e'; validation=run/'results/groups/e2e/03_read_validation'
    qc=read_tsv(run/'results/cohort/isolate_qc.tsv')
    assert len(qc)==2 and all(r['excluded']=='0' for r in qc), qc
    for row in qc:
        assert 0<=float(row['checkm2_completeness'])<=100, row
        assert 0<=float(row['checkm2_contamination'])<=100, row
        verified_archive(validation/'evidence'/row['isolate_id'])
    matrix=validation/'validated_gene_presence_absence.binary.tsv'
    assert read_tsv(matrix), 'Empty validated pangenome'
    assert 'Validation: PASS' in (validation/'summary_statistics.txt').read_text()
    before={str(p):p.stat().st_mtime_ns for p in (run/'state').rglob('*.done.json') if p.parent != run/'state'}
    digest=sha256(matrix)
    subprocess.run([sys.executable,'-m','cleangene','run','--profile','local','--analysis-root',str(root),'--resume','e2e'],check=True)
    assert sha256(matrix)==digest, 'Resume changed completed calls'
    after={str(p):p.stat().st_mtime_ns for p in (run/'state').rglob('*.done.json') if p.parent != run/'state'}
    assert before==after, 'Resume repeated completed stages'
    report={'status':'PASS','run_dir':str(run),'isolates':2,'gene_clusters':len(read_tsv(matrix)),
            'checkm2_evaluated':2,'cram_archives_verified':2,'resume_reused_completed_stages':True,
            'fixture':'200 kb cropped genome; completeness thresholds set to zero for execution testing'}
    (root/'e2e_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
