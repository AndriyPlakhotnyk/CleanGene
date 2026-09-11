"""Downstream read inspection and MSA from retained CRAM evidence."""
from __future__ import annotations

import subprocess
from pathlib import Path

from .alignment_archive import verified_archive, restore_own_bam
from .fasta import read_fasta, write_fasta
from .util import read_tsv, safe_name, write_tsv


def archive_utility(request: dict) -> None:
    from .downstream import load_matrix, subset_samples
    run = Path(request['run_dir']); out = Path(request['output_dir']); group = request['organism']
    root = run / 'results/groups' / safe_name(group) / '03_read_validation'
    isolates, _ = load_matrix(run, group); isolates = subset_samples(isolates, request.get('samples', []))
    out.mkdir(parents=True, exist_ok=True); manifest = []; sequences = {gene: [] for gene in request.get('genes', [])}
    loci = read_tsv(root / 'cluster_isolate_loci.tsv') if (root / 'cluster_isolate_loci.tsv').is_file() else []
    sources = read_tsv(root / 'discovered_gene_sources.tsv') if (root / 'discovered_gene_sources.tsv').is_file() else []
    for isolate in isolates:
        ev = root / 'evidence' / safe_name(isolate); cram, reference, metadata = verified_archive(ev)
        sample = out / safe_name(isolate); sample.mkdir(parents=True, exist_ok=True)
        row = {'isolate_id': isolate, 'cram': str(cram), 'crai': str(cram) + '.crai', 'reference': str(reference), 'reference_sha256': metadata['reference_sha256'], 'bam': '', 'sam': ''}
        kind = request['utility']
        if kind == 'restore_bam': row['bam'] = str(restore_own_bam(ev, sample / 'own_assembly_reads.bam', int(request.get('cpus', 1))))
        elif kind == 'inspect_reads':
            sam = sample / 'reads.sam'
            command = ['samtools', 'view', '-h', '-T', str(reference), str(cram)] + ([request['region']] if request.get('region') else [])
            with sam.open('w') as handle: subprocess.run(command, check=True, stdout=handle)
            row['sam'] = str(sam)
        elif kind == 'evidence_msa':
            from .evidence import consensus, consensus_locus
            # A temporary BAM makes all existing BCFtools options behave identically
            # to validation; it is removed after the consensus has been extracted.
            bam = restore_own_bam(ev, sample / 'temporary.bam', int(request.get('cpus', 1)))
            prefix = sample / 'assembly'
            reconstructed = read_fasta(consensus(reference, bam, prefix, float(request.get('min_depth', 5)), int(request.get('min_mapq', 20)), int(request.get('basequal', 30))))
            for gene in sequences:
                novel = next((r for r in sources if r['Gene'] == gene and r['isolate_id'] == isolate), None)
                if novel:
                    sequence = read_fasta(ev / 'discovered_genes.fasta')[novel['candidate_id']]
                    sequences[gene].append((isolate, sequence)); continue
                candidates = [r for r in loci if r['Gene'] == gene and r['isolate_id'] == isolate]
                for i, locus in enumerate(candidates):
                    sequence = consensus_locus(reconstructed, locus, prefix.with_suffix('.chain'))
                    if sequence: sequences[gene].append((isolate if len(candidates) == 1 else f'{isolate}_copy{i+1}', sequence))
            bam.unlink(); Path(str(bam) + '.bai').unlink(missing_ok=True)
        manifest.append(row)
    write_tsv(out / 'alignment_evidence.tsv', ['isolate_id','cram','crai','reference','reference_sha256','bam','sam'], manifest)
    for gene, records in sequences.items():
        raw = out / f'{safe_name(gene)}.fasta'; aligned = out / f'{safe_name(gene)}.aligned.fasta'
        write_fasta(raw, records)
        if len(records) > 1:
            with aligned.open('w') as handle: subprocess.run(['mafft', '--auto', '--thread', str(request.get('cpus', 1)), str(raw)], check=True, stdout=handle)
        else: write_fasta(aligned, records)
    if sequences: write_tsv(out / 'msa_summary.tsv', ['Gene','sequence_count'], [[gene, len(records)] for gene, records in sequences.items()])
