"""Bounded CDS discovery and deterministic consolidation of validated gene sequences."""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

from .fasta import read_fasta, write_fasta
from .util import read_tsv, run, write_tsv


def arbitration_cases(rows: list[dict], panaroo_present: int, fraction: float = .03) -> tuple[int, list[dict]]:
    if not 0 <= fraction <= 1: raise ValueError('Arbitration fraction must be between 0 and 1')
    maximum = math.floor(panaroo_present * fraction)
    priority = {'divergent_variant': 0, 'possible_truncation': 1, 'partial_homolog': 2, 'insufficient_evidence': 3, 'ambiguous_multimap': 4}
    def number(row, field):
        try: return float(row.get(field) or 0)
        except (ValueError, TypeError): return -1.0
    def key(row):
        state = row.get('evidence_state', '')
        value = number(row, 'identity' if state == 'divergent_variant' else 'breadth')
        return priority.get(state, 5), -value, row['Gene']
    pending = sorted((r for r in rows if r.get('arbitration_status') == 'pending' and r.get('evidence_state') not in {'confirmed_present', 'not_detected'}), key=key)
    return maximum, pending


def complete_cds(seq: str) -> bool:
    from .evidence import orf_integrity
    seq = seq.upper()
    return len(seq) >= 90 and seq[:3] in {'ATG', 'GTG', 'TTG'} and seq[-3:] in {'TAA', 'TAG', 'TGA'} and orf_integrity(seq) == 'intact'


def discover_cds(reconstruction: dict, metric: dict, directory: Path, threads: int,
                 min_depth: float, min_mapq: int, basequal: int) -> dict | None:
    """A full, target-overlapping CDS needs independent recruited-read support.

    Partial homologs retain absence of the original intact gene. A complete CDS
    candidate receives a stable temporary name; consolidation merges equivalences.
    """
    from .evidence import map_reads, coverage, consensus, sequence_identity
    if metric.get('evidence_state') not in {'partial_homolog', 'possible_truncation'}: return None
    match = reconstruction.get('candidate')
    if not match or not reconstruction.get('contigs_path'): return None
    cds = directory / 'predicted_cds.fasta'
    run(['prodigal', '-i', reconstruction['contigs_path'], '-d', str(cds), '-o', str(directory / 'predicted_cds.gff'), '-f', 'gff', '-p', 'meta', '-c', '-q'])
    sequences = read_fasta(cds); candidates = []
    with cds.open() as handle:
        for line in handle:
            if not line.startswith('>') or 'partial=00' not in line: continue
            fields = line[1:].strip().split(' # ')
            if len(fields) < 4: continue
            name = fields[0].split()[0]; contig = name.rsplit('_', 1)[0]
            if contig != match['contig'] or not complete_cds(sequences[name]): continue
            start, end = int(fields[1]) - 1, int(fields[2])
            overlap = max(0, min(end, int(match['query_end'])) - max(start, int(match['query_start'])))
            if overlap >= .5 * max(1, int(match['query_end']) - int(match['query_start'])):
                candidates.append((overlap, name, sequences[name]))
    if not candidates: return None
    _, _, sequence = max(candidates, key=lambda x: (x[0], len(x[2]), x[1]))
    gene = 'CGNEW_' + hashlib.sha256(sequence.encode()).hexdigest()[:20]
    reference = directory / 'novel_candidate.fasta'; write_fasta(reference, [(gene, sequence)])
    run(['bwa', 'index', str(reference)], stdout=directory / 'novel_index.stdout', stderr=directory / 'novel_index.stderr')
    run(['samtools', 'faidx', str(reference)])
    # Recruitment retained both mates. Also include singleton/other reads.
    read1, read2 = directory / 'novel_R1.fastq', directory / 'novel_R2.fastq'
    import shutil
    with read1.open('wb') as output:
        for source in (directory / 'recruited_R1.fastq', directory / 'all_singletons.fastq'):
            if source.is_file():
                with source.open('rb') as handle: shutil.copyfileobj(handle, output)
    shutil.copyfile(directory / 'recruited_R2.fastq', read2)
    # Interleaving singleton records in R1 would unbalance pairing; use single-end
    # validation of all recruited reads for this candidate, preserving every base.
    combined = directory / 'novel_reads.fastq'
    with combined.open('wb') as output:
        for source in (read1, read2):
            with source.open('rb') as handle: shutil.copyfileobj(handle, output)
    bam = directory / 'novel_reads.bam'
    map_reads(reference, str(combined), '', bam, threads, min_mapq, directory / 'novel_bwa.log', retain_ambiguous=True)
    cov = coverage(bam, min_mapq).get(gene, {})
    reconstructed = read_fasta(consensus(reference, bam, directory / 'novel', min_depth, min_mapq, basequal)).get(gene, '')
    identity = sequence_identity(sequence, reconstructed)
    if not (cov.get('breadth', 0) >= .95 and cov.get('mean_depth', 0) >= min_depth and identity and identity['identity'] >= .95 and complete_cds(reconstructed)):
        return None
    gene = 'CGNEW_' + hashlib.sha256(reconstructed.encode()).hexdigest()[:20]
    return {'candidate_id': gene, 'parent_gene': metric['Gene'], 'sequence': reconstructed,
            'evidence_state': 'discovered_complete_cds', 'breadth': cov['breadth'], 'identity': identity['identity'],
            'mean_depth': cov['mean_depth'], 'mapped_reads': cov['mapped_reads'], 'orf_integrity': 'intact',
            'discovery_reason': metric['evidence_state'], 'sequence_sha256': hashlib.sha256(reconstructed.encode()).hexdigest()}


def cluster_new_sequences(records: list[tuple[str, str]], directory: Path, identity: float = .95, coverage: float = .95, threads: int = 1) -> tuple[dict[str, str], dict[str, str]]:
    """CD-HIT-EST global nucleotide identity with reciprocal length coverage."""
    if not .8 <= identity <= 1 or not 0 < coverage <= 1: raise ValueError('Invalid nucleotide clustering thresholds')
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / 'new_genes.fasta'; output = directory / 'new_genes.clustered.fasta'
    write_fasta(fasta, sorted(records, key=lambda x: (-len(x[1]), x[0])))
    if not records: write_fasta(output, []); return {}, {}
    run(['cd-hit-est', '-i', str(fasta), '-o', str(output), '-c', str(identity), '-G', '1', '-aS', str(coverage), '-aL', str(coverage), '-n', '8', '-d', '0', '-T', str(threads), '-M', '0'], stdout=directory / 'clustering.stdout', stderr=directory / 'clustering.stderr')
    clusters = []; members = []; representative = ''
    for line in Path(str(output) + '.clstr').read_text().splitlines():
        if line.startswith('>Cluster'):
            if members: clusters.append((representative, members))
            members = []; representative = ''
        else:
            name = line.split('>', 1)[1].split('...', 1)[0]; members.append(name)
            if line.rstrip().endswith('*'): representative = name
    if members: clusters.append((representative, members))
    sequences = dict(records); aliases = {}; representatives = {}
    for representative, members in clusters:
        if not representative: raise ValueError('CD-HIT cluster lacks representative')
        sequence = sequences[representative]
        name = 'CGNEW_' + hashlib.sha256(sequence.encode()).hexdigest()[:20]
        representatives[name] = sequence
        for member in members: aliases[member] = name
    if set(aliases) != set(sequences): raise ValueError('Incomplete novel-gene clustering')
    return representatives, aliases


def consolidate_discoveries(out: Path, isolates: list[str], existing: dict[str, str], threads: int = 1,
                           identity: float = .95, coverage: float = .95) -> tuple[dict[str,str], dict[str,dict[str,int]], list[dict]]:
    """Consolidate discoveries after every isolate has finished arbitration.

    Existing Panaroo names remain stable. Only supported discoveries receive 1;
    other new-gene cells are explicitly untested, not confirmed absence.
    """
    discoveries = []; seqs = {}
    from .util import safe_name
    for isolate in isolates:
        directory = out / 'evidence' / safe_name(isolate)
        fasta = directory / 'discovered_genes.fasta'; table = directory / 'discovered_genes.tsv'
        if not table.is_file(): continue
        sequences = read_fasta(fasta)
        for row in read_tsv(table):
            seqs[row['candidate_id']] = sequences[row['candidate_id']]
            discoveries.append({**row, 'isolate_id': isolate})
    representatives, aliases = cluster_new_sequences(list(seqs.items()), out / 'consolidation', identity, coverage, threads)
    # Compare discovered representatives to the existing reference catalogue. Full
    # reciprocal coverage avoids merging a short homolog into a longer gene.
    forbidden = {(aliases[row["candidate_id"]], row["parent_gene"]) for row in discoveries
                 if row.get("discovery_reason") == "partial_homolog"
                 and float(row.get("parent_read_breadth") or 0) < float(row.get("parent_min_breadth") or .95)}
    rename = {}
    if representatives and existing:
        directory = out / 'consolidation'
        ref = directory / 'existing.fasta'; query = directory / 'representatives.fasta'
        existing_ids = {f'EXISTING_{i}': gene for i, gene in enumerate(sorted(existing))}
        write_fasta(ref, [(name, existing[gene]) for name, gene in existing_ids.items()])
        write_fasta(query, list(representatives.items()))
        output = directory / 'distinct_new_genes.fasta'
        run(['cd-hit-est-2d', '-i', str(ref), '-i2', str(query), '-o', str(output),
             '-c', str(identity), '-G', '1', '-aS', str(coverage), '-aL', str(coverage),
             '-s2', str(coverage), '-S2', '999999', '-n', '8', '-d', '0', '-g', '1',
             '-T', str(threads), '-M', '0'], stdout=directory / 'consolidation.stdout', stderr=directory / 'consolidation.stderr')
        cluster = []
        def register(members):
            reference = next((existing_ids[name] for name in members if name in existing_ids), None)
            if reference:
                for name in members:
                    if name in representatives and (name, reference) not in forbidden: rename[name] = reference
        for line in Path(str(output) + '.clstr').read_text().splitlines():
            if line.startswith('>Cluster'):
                register(cluster); cluster = []
            else: cluster.append(line.split('>', 1)[1].split('...', 1)[0])
        register(cluster)
    merged = {name: seq for name, seq in representatives.items() if name not in rename}
    calls = {}; sources = []
    for row in discoveries:
        representative = aliases[row['candidate_id']]; gene = rename.get(representative, representative)
        calls.setdefault(gene, {iso: 0 for iso in isolates})[row['isolate_id']] = 1
        sources.append({**row, 'Gene': gene, 'representative_id': representative})
    fields = ['Gene','isolate_id','candidate_id','parent_gene','representative_id','discovery_reason','sequence_sha256','breadth','identity','mean_depth','mapped_reads','orf_integrity']
    write_tsv(out / 'discovered_gene_sources.tsv', fields, sources)
    write_tsv(out / 'discovered_gene_aliases.tsv', ['candidate_id','Gene'], [[member, rename.get(rep, rep)] for member, rep in sorted(aliases.items())])
    write_fasta(out / 'validated_gene_sequences.fasta', [*existing.items(), *merged.items()])
    return merged, calls, sources
