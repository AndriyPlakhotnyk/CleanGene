import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cleangene.alignment_archive import (alignment_digest, archive_own_alignment,
    retain_reference, restore_own_bam, sha256, verified_archive)
from cleangene.discovery import arbitration_cases, complete_cds, consolidate_discoveries, discover_cds
from cleangene.evidence import classify_gene_evidence
from cleangene.fasta import read_fasta, write_fasta
from cleangene.util import read_tsv, write_tsv


def have(*tools):
    return all(shutil.which(tool) for tool in tools)


def archive_fixture(root):
    ev = root / 'evidence'; ev.mkdir(parents=True)
    source = root / 'assembly.fa'; source.write_text('>ctg exact reference\n' + 'ACGT' * 100 + '\n')
    reference = retain_reference(source, ev)
    sam = root / 'reads.sam'
    sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ctg\tLN:400\n'
        'r1\t0\tctg\t1\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\tNM:i:0\tMD:Z:20\tZZ:Z:custom\n'
        'r2\t256\tctg\t5\t0\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tHHHHHHHHHHHHHHHHHHHH\tAS:i:20\n'
        'r3\t4\t*\t0\t0\t*\t*\t0\t0\tNNACGT\t!!IIII\tXY:B:i,1,2,3\n')
    bam = ev / 'own_assembly_reads.bam'
    subprocess.run(['samtools', 'view', '-b', '-o', str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run(['samtools', 'index', str(bam)], check=True, capture_output=True)
    return ev, reference, bam


class DecisionRoutingTests(unittest.TestCase):
    def test_terminal_calls_and_depth_gate(self):
        for initial in (0, 1):
            zero = classify_gene_evidence(initial_call=initial, mapped_reads=0, breadth=0, mean_depth=0, identity=None)
            self.assertEqual((zero['validated_call'], zero['final_call_source']), (0, 'read_validation'))
            strong = classify_gene_evidence(initial_call=initial, mapped_reads=50, breadth=.95, mean_depth=5, identity=.95)
            self.assertEqual(strong['validated_call'], 1)
            self.assertNotEqual(strong['final_call_source'], 'arbitration_pending')
            low = classify_gene_evidence(initial_call=initial, mapped_reads=50, breadth=1, mean_depth=4.99, identity=1)
            self.assertEqual(low['final_call_source'], 'arbitration_pending')
        ambiguous = classify_gene_evidence(mapped_reads=10, breadth=0, mean_depth=0, identity=None)
        self.assertNotEqual(ambiguous['evidence_state'], 'not_detected')

    def test_three_percent_and_descending_state_ranking(self):
        rows = [{'Gene': name, 'evidence_state': state, 'arbitration_status': 'pending', 'identity': identity, 'breadth': breadth}
                for name, state, identity, breadth in [
                    ('p_low', 'partial_homolog', 1, .3), ('d_low', 'divergent_variant', .90, 1),
                    ('t_high', 'possible_truncation', .97, .94), ('d_high', 'divergent_variant', .94, .91),
                    ('t_low', 'possible_truncation', 1, .8), ('p_high', 'partial_homolog', .8, .9),
                    ('i', 'insufficient_evidence', '', 1), ('zero', 'not_detected', '', 0),
                    ('strong', 'confirmed_present', 1, 1)]]
        cap, cases = arbitration_cases(rows, 100)
        self.assertEqual(cap, 3)
        self.assertEqual([r['Gene'] for r in cases], ['d_high','d_low','t_high','t_low','p_high','p_low','i'])
        self.assertEqual(arbitration_cases(rows, 33)[0], 0)
        self.assertEqual(arbitration_cases(rows, 34)[0], 1)
        self.assertEqual(arbitration_cases(rows, 0)[0], 0)


@unittest.skipUnless(have('samtools'), 'SAMtools required')
class ArchiveTests(unittest.TestCase):
    def test_lossless_roundtrip_reference_and_idempotency(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); ev, reference, bam = archive_fixture(root)
            expected = alignment_digest(bam, reference)
            metadata = archive_own_alignment(ev)
            self.assertFalse(bam.exists())
            self.assertFalse(Path(str(bam) + '.bai').exists())
            self.assertTrue((ev / 'own_assembly_reads.cram.crai').is_file())
            self.assertEqual(metadata['reference_sha256'], sha256(root / 'assembly.fa'))
            self.assertEqual(metadata['alignment_count'], 3)
            self.assertEqual(archive_own_alignment(ev), metadata)
            restored = restore_own_bam(ev, root / 'restored.bam')
            self.assertEqual(alignment_digest(restored, reference), expected)
            self.assertEqual(restore_own_bam(ev, restored), restored)
            self.assertTrue((ev / 'own_assembly_reads.cram').is_file())
            reference.write_text('>ctg\n' + 'TGCA' * 100 + '\n')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'): verified_archive(ev)

    def test_conversion_failure_keeps_original_bam(self):
        with tempfile.TemporaryDirectory() as d:
            ev, reference, bam = archive_fixture(Path(d))
            before = bam.read_bytes()
            with patch('cleangene.alignment_archive.run', side_effect=subprocess.CalledProcessError(1, 'samtools')):
                with self.assertRaises(subprocess.CalledProcessError): archive_own_alignment(ev)
            self.assertEqual(bam.read_bytes(), before)
            self.assertFalse((ev / 'own_assembly_reads.archive.json').exists())

    def test_changed_reference_cannot_replace_archive_reference(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); ev, reference, bam = archive_fixture(root)
            archive_own_alignment(ev); before = reference.read_bytes()
            source = root / 'other.fa'; source.write_text('>other\nAAAA\n')
            with self.assertRaisesRegex(ValueError, 'Assembly differs'): retain_reference(source, ev)
            self.assertEqual(reference.read_bytes(), before)
            verified_archive(ev)

    def test_read_inspection_and_restoration_utilities(self):
        from cleangene.archive_utils import archive_utility
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); group = root / 'results/groups/g/03_read_validation'
            ev, reference, bam = archive_fixture(group)
            ev.rename(group / 'tmp_evidence'); (group / 'evidence').mkdir()
            (group / 'tmp_evidence').rename(group / 'evidence/i')
            ev = group / 'evidence/i'; archive_own_alignment(ev)
            request = {'run_dir': str(root), 'organism': 'g', 'samples': [], 'output_dir': str(root / 'inspect'), 'utility': 'inspect_reads', 'region': 'ctg:1-30'}
            with patch('cleangene.downstream.load_matrix', return_value=(['i'], {})):
                archive_utility(request)
                self.assertIn('r1\t', (root / 'inspect/i/reads.sam').read_text())
                self.assertNotIn('r3\t', (root / 'inspect/i/reads.sam').read_text())
                archive_utility({**request, 'utility': 'restore_bam', 'output_dir': str(root / 'restore')})
            self.assertTrue((root / 'restore/i/own_assembly_reads.bam.bai').exists())


@unittest.skipUnless(have('cd-hit-est', 'cd-hit-est-2d'), 'CD-HIT required')
class ConsolidationTests(unittest.TestCase):
    def test_cluster_discoveries_merge_existing_and_preserve_distinct(self):
        rng = random.Random(42)
        original = ''.join(rng.choice('ACGT') for _ in range(1200))
        variant = ''.join(('A' if base != 'A' else 'C') if i % 50 == 0 else base for i, base in enumerate(original))
        known = ''.join(rng.choice('ACGT') for _ in range(120))
        known_variant = ('A' if known[0] != 'A' else 'C') + known[1:]
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for iso, records in [('a', [('candidateA', original), ('rediscovered', known_variant)]), ('b', [('candidateB', variant)])]:
                ev = root / 'evidence' / iso
                write_fasta(ev / 'discovered_genes.fasta', records)
                write_tsv(ev / 'discovered_genes.tsv', ['candidate_id','parent_gene','discovery_reason'], [[name, 'parent', 'partial_homolog'] for name, seq in records])
            merged, calls, sources = consolidate_discoveries(root, ['a','b','c'], {'known': known})
            self.assertEqual(len(merged), 1)
            gene = next(iter(merged))
            self.assertEqual(calls[gene], {'a':1,'b':1,'c':0})
            self.assertEqual(calls['known'], {'a':1,'b':0,'c':0})
            self.assertEqual(set(read_fasta(root / 'validated_gene_sequences.fasta')), {'known',gene})
            self.assertEqual(len(sources), 3)
            self.assertEqual(len(read_tsv(root / 'discovered_gene_aliases.tsv')), 3)


@unittest.skipUnless(have('prodigal','bwa','samtools','bcftools'), 'CDS discovery tools required')
class CDSDiscoveryTests(unittest.TestCase):
    def test_complete_cds_requires_read_supported_full_orf(self):
        rng = random.Random(81)
        codons = ['GCT','GCC','GAA','GAC','AAA','AAG','CTG','ATC','GGT','ACC','TTC','CAG','CGT','GTG']
        sequence = 'ATG' + ''.join(rng.choice(codons) for _ in range(500)) + 'TAA'
        flank = ''.join(rng.choice('ACGT') for _ in range(300))
        contig = flank + sequence + flank
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); assembly = root / 'contigs.fasta'; write_fasta(assembly, [('ctg',contig)])
            with (root / 'recruited_R1.fastq').open('w') as handle:
                for i, start in enumerate(range(0, len(contig)-150, 3)):
                    handle.write(f'@r{i}\n{contig[start:start+150]}\n+\n' + 'I'*150 + '\n')
            (root / 'recruited_R2.fastq').touch()
            reconstruction = {'contigs_path':str(assembly), 'candidate':{'contig':'ctg','query_start':600,'query_end':1200}}
            metric = {'Gene':'parent','evidence_state':'partial_homolog'}
            discovered = discover_cds(reconstruction, metric, root, 1, 5, 20, 30)
            self.assertIsNotNone(discovered)
            self.assertTrue(complete_cds(discovered['sequence']))
            self.assertEqual(discovered['parent_gene'], 'parent')
            self.assertGreaterEqual(discovered['breadth'], .95)
            self.assertGreaterEqual(discovered['identity'], .95)
            self.assertTrue(discovered['candidate_id'].startswith('CGNEW_'))
            # A predicted complete CDS alone is insufficient without read depth.
            with patch('cleangene.evidence.coverage', return_value={}):
                self.assertIsNone(discover_cds(reconstruction, metric, root, 1, 5, 20, 30))


@unittest.skipUnless(have('samtools','bcftools','mafft'), 'MSA tools required')
class ArchivedMSATests(unittest.TestCase):
    def test_msa_reads_archives_and_removes_temporary_bams(self):
        from cleangene.archive_utils import archive_utility
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); group = root / 'results/groups/g/03_read_validation'
            for isolate in ('a','b'):
                ev, ref, bam = archive_fixture(root / isolate)
                destination = group / 'evidence' / isolate
                destination.parent.mkdir(parents=True, exist_ok=True); ev.rename(destination)
                archive_own_alignment(destination)
            write_tsv(group / 'cluster_isolate_loci.tsv', ['Gene','isolate_id','assembly_scaffold','cds_start','cds_end','cds_strand'],
                      [['gene', iso, 'ctg', 1, 20, '+'] for iso in ('a','b')])
            out = root / 'msa'
            request = {'run_dir':str(root),'output_dir':str(out),'organism':'g','utility':'evidence_msa','genes':['gene'],'min_depth':1,'cpus':1}
            with patch('cleangene.downstream.load_matrix', return_value=(['a','b'], {})): archive_utility(request)
            self.assertEqual(read_fasta(out / 'gene.aligned.fasta'), {'a':'ACGT'*5,'b':'ACGT'*5})
            self.assertEqual(read_tsv(out / 'msa_summary.tsv')[0]['sequence_count'], '2')
            self.assertFalse((out / 'a/temporary.bam').exists())
            self.assertTrue((group / 'evidence/a/own_assembly_reads.cram').exists())


@unittest.skipUnless(have('cd-hit-est','cd-hit-est-2d'), 'CD-HIT required')
class ReductionDiscoveryTests(unittest.TestCase):
    def test_reduction_adds_cluster_calls_and_validated_summary(self):
        from cleangene.util import atomic_json
        from cleangene.workers import reduce_group
        rng = random.Random(9); sequence = ''.join(rng.choice('ACGT') for _ in range(900))
        with tempfile.TemporaryDirectory() as d:
            run = Path(d); root = run / 'results/groups/g'; out = root / '03_read_validation'
            atomic_json(run / 'provenance/resolved_config.json', {})
            write_tsv(root / '02_pangenome/initial_calls/gene_presence_absence.binary.tsv', ['Gene','a','b'], [['parent',1,0]])
            write_fasta(out / 'tested_gene_references.fasta', [('ref','A'*900)])
            write_tsv(out / 'tested_gene_key.tsv', ['Gene','reference_id'], [['parent','ref']])
            write_fasta(out / 'evidence/a/discovered_genes.fasta', [('new',sequence)])
            write_tsv(out / 'evidence/a/discovered_genes.tsv', ['candidate_id','parent_gene','discovery_reason'], [['new','parent','partial_homolog']])
            write_tsv(out / 'evidence/a/arbitrated_metrics.tsv', ['Gene','validated_call','evidence_state'], [['parent',0,'partial_homolog']])
            with patch('cleangene.workers.task_row', return_value={'group_id':'g'}), patch('cleangene.workers.retained_rows', return_value=[{'isolate_id':'a'},{'isolate_id':'b'}]):
                reduce_group(run,0)
            matrix = read_tsv(out / 'validated_gene_presence_absence.binary.tsv')
            self.assertEqual(len(matrix),2)
            self.assertEqual((matrix[1]['a'],matrix[1]['b']),('1','0'))
            self.assertEqual(matrix[0]['a'],'0')
            evidence = read_tsv(out / 'gene_call_evidence.long.tsv')
            self.assertEqual(evidence[-1]['evidence_state'],'new_gene_not_tested')
            summary = read_tsv(out / 'gene_call_summary.tsv')[0]
            self.assertEqual((summary['total_calls'],summary['changed_0_to_1'],summary['changed_1_to_0']),('4','1','1'))
            self.assertIn('Validation: PASS',(out / 'summary_statistics.txt').read_text())


class ArchiveCLIWorkflowTests(unittest.TestCase):
    def test_local_execution_and_slurm_submission_share_request(self):
        from cleangene.utils_cli import main
        from cleangene.util import atomic_json
        with tempfile.TemporaryDirectory() as d:
            run = Path(d)
            atomic_json(run / 'provenance/resolved_config.json', {'UTILS_CPUS':'1'})
            write_tsv(run / 'provenance/manifest.tsv', ['isolate_id','group_id'], [['i','g']])
            write_tsv(run / 'state/group_tasks.tsv', ['group_id'], [['g']])
            write_tsv(run / 'results/groups/g/cleaned_pangenome.tsv', ['Gene','i'], [['gene',1]])
            with patch('cleangene.downstream.run_request') as execute, patch('cleangene.utils_cli.submit',return_value='123') as submit:
                self.assertEqual(main(['restore-bam','--run-dir',str(run),'--organism','g','--profile','local','--analysis-name','local']),0)
                execute.assert_called_once(); submit.assert_not_called()
                self.assertEqual(main(['evidence-msa','--run-dir',str(run),'--organism','g','--genes','gene','--profile','slurm','--analysis-name','slurm']),0)
                submit.assert_called_once()
                self.assertIn('_utils_worker', ' '.join(submit.call_args.args[0]))


@unittest.skipUnless(have('samtools'), 'SAMtools required')
class ArbitrationWorkerTests(unittest.TestCase):
    def test_cap_order_terminal_skip_and_archive_before_completion(self):
        from cleangene.util import atomic_json, load_json
        from cleangene.workers import arbitrate
        with tempfile.TemporaryDirectory() as d:
            run = Path(d); group = run / 'results/groups/g'; out = group / '03_read_validation'
            ev, reference, bam = archive_fixture(out)
            ev.rename(out / 'staging'); (out / 'evidence').mkdir(); (out / 'staging').rename(out / 'evidence/i')
            ev = out / 'evidence/i'
            atomic_json(run / 'provenance/resolved_config.json', {})
            write_tsv(group / '02_pangenome/initial_calls/gene_presence_absence.binary.tsv', ['Gene','i'], [[f'g{n}',1] for n in range(100)])
            fields = ['Gene','initial_call','validated_call','evidence_state','arbitration_status','identity','breadth','reference_id']
            write_tsv(ev / 'metrics.tsv', fields, [
                ['low',1,0,'partial_homolog','pending',.99,.5,'low'],
                ['d_low',1,1,'divergent_variant','pending',.91,.95,'d_low'],
                ['t',1,1,'possible_truncation','pending',.99,.94,'t'],
                ['d_high',1,1,'divergent_variant','pending',.94,.95,'d_high'],
                ['zero',1,0,'not_detected','not_required','',0,'zero'],
                ['strong',0,1,'confirmed_present','not_required',1,1,'strong'],
            ])
            with patch('cleangene.workers.task_row',return_value={'group_id':'g','isolate_id':'i'}), patch('cleangene.workers.retained_rows',return_value=[{'isolate_id':'i'}]), patch('cleangene.workers.targeted_local_reconstruction',return_value={'status':'no_recruited_reads'}) as reconstruct:
                arbitrate(run,0)
            self.assertEqual([c.kwargs['region'] for c in reconstruct.call_args_list],['d_high','d_low','t'])
            metadata = load_json(run / 'state/arbitrate/i.done.json')
            self.assertEqual((metadata['arbitration_cap'],metadata['processed_cases'],metadata['deferred_cases']),(3,3,1))
            rows = {r['Gene']:r for r in read_tsv(ev / 'arbitrated_metrics.tsv')}
            self.assertEqual(rows['low']['arbitration_status'],'deferred_limit')
            self.assertEqual(rows['zero']['validated_call'],'0')
            self.assertEqual(rows['strong']['validated_call'],'1')
            self.assertFalse((ev / 'own_assembly_reads.bam').exists())
            verified_archive(ev)
            with patch('cleangene.workers.task_row',return_value={'group_id':'g','isolate_id':'i'}), patch('cleangene.workers.targeted_local_reconstruction') as reconstruct:
                arbitrate(run,0)
                reconstruct.assert_not_called()
