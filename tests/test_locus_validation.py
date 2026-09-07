import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from cleangene.evidence import (EVIDENCE_VERSION, classify_gene_evidence, consensus_locus,
    coverage, map_reads, orf_integrity, sequence_identity, spans_junction,
    targeted_local_reconstruction, validate_isolate)
from cleangene.util import read_tsv, write_tsv
from cleangene.workers import arbitrate_evidence


class LocusValidationTests(unittest.TestCase):
    def test_high_identity_93_percent_is_possible_truncation(self):
        result=classify_gene_evidence(initial_call=1,mapped_reads=50,breadth=.93,mean_depth=20,identity=.99)
        self.assertEqual((result['evidence_state'],result['validated_call']),('possible_truncation',1))

    def test_masked_orf_is_unresolved(self):
        self.assertEqual(orf_integrity('ATGNNNTAA'),'unresolved')

    def test_indel_comparison_does_not_shift_downstream_identity(self):
        result=sequence_identity('ACGTACGTACGT','ACGTATCGTACGT')
        self.assertAlmostEqual(result['identity'],12/13)
        self.assertIsNone(sequence_identity('ATG','NNN'))

    def test_lifts_cds_after_insertion_on_both_strands(self):
        with tempfile.TemporaryDirectory() as d:
            chain=Path(d)/'own.chain'; chain.write_text('chain 0 ctg 9 + 0 9 ctg 11 + 0 11 1\n3 0 2\n6\n')
            row={'assembly_scaffold':'ctg','cds_start':'4','cds_end':'9','cds_strand':'+'}
            self.assertEqual(consensus_locus({'ctg':'AAATTCCCGGG'},row,chain),'TTCCCGGG')
            row['cds_start']='7'; row['cds_strand']='-'
            self.assertEqual(consensus_locus({'ctg':'AAATTCCCGGG'},row,chain),'CCC')

    def test_deletion_requires_contiguous_both_flank_anchors(self):
        self.assertTrue(spans_junction({'reference_start':0,'alignment_cigar':'1000M'},500))
        self.assertFalse(spans_junction({'reference_start':0,'alignment_cigar':'500M200I500M'},500))
        self.assertFalse(spans_junction({'reference_start':0,'alignment_cigar':'450M100D450M'},500))
        self.assertFalse(spans_junction({'reference_start':0,'alignment_cigar':'500M'},500))

    def test_family_call_stays_unresolved_after_arbitration(self):
        for initial in ('0','1'):
            result=arbitrate_evidence({'initial_call':initial,'evidence_state':'ambiguous_multimap'})
            self.assertEqual(result['validated_call'],'')
            self.assertEqual(result['arbitration_status'],'family_only')

    def test_secondary_only_coverage_is_retained(self):
        with patch('cleangene.evidence.subprocess.run',return_value=subprocess.CompletedProcess([],0,'g\t1\t100\t20\t100\t100\t20\t30\t0\n')) as run:
            result=coverage(Path('reads.bam'),20,include_ambiguous=True)
            self.assertEqual(result['g']['mapped_reads'],20)
            self.assertIn('1540',run.call_args.args[0])

    def test_matching_mapping_marker_reuses_bam(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); paths=[root/x for x in ('ref.fa','r1.fq','r2.fq')]
            for p in paths: p.write_text('data')
            bam=root/'reads.bam'; bam.touch(); Path(str(bam)+'.bai').touch()
            signature={'version':EVIDENCE_VERSION,'min_mapq':20,'retain_ambiguous':True,'inputs':[[str(p.resolve()),p.stat().st_size,p.stat().st_mtime_ns] for p in paths]}
            bam.with_suffix('.mapping.json').write_text(json.dumps(signature))
            with patch('cleangene.evidence.subprocess.Popen') as process:
                map_reads(*paths,bam,2,20,root/'map.log',retain_ambiguous=True)
                process.assert_not_called()

    def test_reconstruction_collates_mates_and_creates_match_directory(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); commands=[]
            def run(cmd,**kwargs):
                commands.append(cmd)
                if cmd[:2]==['samtools','fastq']:
                    for flag in ('-1','-2','-0','-s'): Path(cmd[cmd.index(flag)+1]).write_text('@r\nACGT\n+\nIIII\n')
                if cmd[0]=='spades.py':
                    out=Path(cmd[cmd.index('-o')+1]); out.mkdir(); (out/'contigs.fasta').write_text('>ctg\nACGT\n')
            def subprocess_run(cmd,**kwargs):
                if cmd[0]=='minimap2':
                    self.assertTrue(Path(cmd[-2]).is_file())
                    return subprocess.CompletedProcess(cmd,0,'ctg\t4\t0\t4\t+\ttarget\t4\t0\t4\t4\t4\t60\tcg:Z:4M\n')
                return subprocess.CompletedProcess(cmd,0,'r\t0\n')
            import io
            from unittest.mock import MagicMock
            process=MagicMock(); process.stdout=io.StringIO('r\t0\n'); process.wait.return_value=0; process.__enter__.return_value=process
            with patch('cleangene.evidence.run',side_effect=run), patch('cleangene.evidence.subprocess.run',side_effect=subprocess_run), patch('cleangene.evidence.subprocess.Popen',return_value=process):
                result=targeted_local_reconstruction(bam=root/'reads.bam',region='g',reference_seq='ACGT',outdir=root,threads=2)
            self.assertEqual(result['status'],'reconstructed')
            fastq=next(c for c in commands if c[:2]==['samtools','fastq'])
            self.assertTrue(fastq[-1].endswith('recruited.collated.bam'))
            self.assertEqual((root/'all_singletons.fastq').read_text().count('@r'),2)

    def test_own_locus_positive_ignores_competitive_zero(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); assembly=root/'assembly.fa'; assembly.write_text('>ctg\nATGAAATAA\n')
            ref=root/'ref.fa'; ref.write_text('>r\nATGAAATAA\n')
            key=root/'key.tsv'; write_tsv(key,['reference_id','Gene'],[['r','g']])
            loci=root/'loci.tsv'; write_tsv(loci,['Gene','assembly_scaffold','cds_start','cds_end','cds_strand'],[['g','ctg',1,9,'+']])
            def consensus(ref,bam,prefix,*args):
                prefix.with_suffix('.chain').write_text('chain 0 ctg 9 + 0 9 ctg 9 + 0 9 1\n9\n')
                return ref
            metric={('ctg',1,9):{'mapped_reads':50,'ambiguous_mapped_reads':0,'breadth':1.,'mean_depth':20.}}
            with patch('cleangene.evidence.run'),patch('cleangene.evidence.map_reads') as mapper,patch('cleangene.evidence.coverage',return_value={}),patch('cleangene.evidence.consensus',side_effect=consensus),patch('cleangene.evidence.locus_coverage',return_value=metric):
                validate_isolate(ref,key,loci,assembly,'r1','r2',root,2,.95,5,.95,20,30,initial_calls={'g':1})
            row=read_tsv(root/'metrics.tsv')[0]
            self.assertEqual(row['validated_call'],'1')
            self.assertEqual(row['final_call_source'],'own_locus_read_validation')
            self.assertEqual(mapper.call_count,2)

    def test_local_dispatches_same_arbitration_worker(self):
        from cleangene.cli import local
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); (root/'state').mkdir(); (root/'provenance').mkdir()
            (root/'provenance'/'resolved_config.json').write_text('{}')
            write_tsv(root/'provenance'/'manifest.tsv',['isolate_id'],[['i']])
            write_tsv(root/'state'/'isolate_tasks.tsv',['isolate_id'],[['i']])
            write_tsv(root/'state'/'group_tasks.tsv',['group_id'],[['g']])
            with patch('cleangene.cli.needs_kraken',return_value=False),patch('cleangene.cli.needs_checkm2',return_value=False),patch('cleangene.cli.dispatch') as dispatch:
                local(root)
            stages=[call.args[0] for call in dispatch.call_args_list]
            self.assertLess(stages.index('validate'),stages.index('arbitrate'))
            self.assertLess(stages.index('arbitrate'),stages.index('reduce'))

    def test_partial_homolog_does_not_become_present_without_reconstruction(self):
        row={'initial_call':'1','validated_call':'0','evidence_state':'partial_homolog','arbitration_status':'pending'}
        result=arbitrate_evidence(row)
        self.assertEqual(result['validated_call'],'0')
        self.assertEqual(result['evidence_state'],'partial_homolog')

    def test_two_bam_passes_measure_overlapping_cds_and_secondary_reads(self):
        import io
        from unittest.mock import MagicMock
        from cleangene.evidence import locus_coverage
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            loci=[{'assembly_scaffold':'ctg','cds_start':'1','cds_end':'3'}, {'assembly_scaffold':'ctg','cds_start':'3','cds_end':'5'}]
            outputs=['ctg\t1\t2\nctg\t3\t4\n', 'r\t0\tctg\t1\t60\t3M\nq\t256\tctg\t3\t0\t3M\n']
            def process(*args,**kwargs):
                p=MagicMock(); p.stdout=io.StringIO(outputs.pop(0)); p.wait.return_value=0; p.__enter__.return_value=p; return p
            with patch('cleangene.evidence.subprocess.Popen',side_effect=process) as popen:
                result=locus_coverage(root/'bam',loci,root,20)
            self.assertEqual(popen.call_count,2)
            self.assertAlmostEqual(result[('ctg',1,3)]['breadth'],2/3)
            self.assertEqual(result[('ctg',3,5)]['ambiguous_mapped_reads'],1)
            self.assertEqual(result[('ctg',3,5)]['mapped_reads'],1)

    def test_resume_invalidates_old_evidence_but_preserves_new_metrics(self):
        from cleangene.workers import invalidate_legacy_identity_metrics
        from cleangene.defaults import DEFAULTS
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); metrics=root/'results/groups/g/03_read_validation/evidence/i/metrics.tsv'
            write_tsv(metrics,['evidence_version','Gene'],[['1','g']])
            self.assertEqual(invalidate_legacy_identity_metrics(root,DEFAULTS),1)
            write_tsv(metrics,['evidence_version','Gene'],[[EVIDENCE_VERSION,'g']])
            self.assertEqual(invalidate_legacy_identity_metrics(root,DEFAULTS),0)

    @unittest.skipUnless(all(__import__('shutil').which(tool) for tool in ('bwa','samtools','bcftools','minimap2','spades.py')), 'BWA/SAMtools/BCFtools are not on PATH')
    def test_real_reads_own_locus_ambiguous_homolog_and_missed_annotation(self):
        import random
        random_source=random.Random(17)
        sequence=''.join(random_source.choice('ACGT') for _ in range(4000))
        complement=str.maketrans('ACGT','TGCA')
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); assembly=root/'assembly.fa'; assembly.write_text('>ctg\n'+sequence+'\n')
            ref=root/'genes.fa'; ref.write_text('>r1\n'+sequence[600:1500]+'\n>r2\n'+sequence[600:1500]+'\n>r3\n'+sequence[2300:3200]+'\n')
            r1=root/'r1.fq'; r2=root/'r2.fq'
            with r1.open('w') as left,r2.open('w') as right:
                for index,start in enumerate(range(0,3600,5)):
                    left.write(f'@read{index}\n{sequence[start:start+150]}\n+\n'+('I'*150)+'\n')
                    right.write(f'@read{index}\n{sequence[start+250:start+400].translate(complement)[::-1]}\n+\n'+('I'*150)+'\n')
            key=root/'key.tsv'; write_tsv(key,['reference_id','Gene'],[['r1','present'],['r2','homolog'],['r3','missed']])
            loci=root/'loci.tsv'; write_tsv(loci,['Gene','assembly_scaffold','cds_start','cds_end','cds_strand'],[['present','ctg',601,1500,'+']])
            subprocess.run(['bwa','index',str(ref)],check=True,capture_output=True)
            subprocess.run(['samtools','faidx',str(ref)],check=True,capture_output=True)
            validate_isolate(ref,key,loci,assembly,str(r1),str(r2),root/'evidence',2,.95,5,.95,20,30,initial_calls={'present':1,'homolog':0,'missed':0})
            rows={row['Gene']:row for row in read_tsv(root/'evidence/metrics.tsv')}
            self.assertEqual(rows['present']['evidence_state'],'confirmed_present')
            self.assertEqual(rows['homolog']['evidence_state'],'ambiguous_multimap')
            self.assertEqual(rows['homolog']['validated_call'],'')
            self.assertEqual(rows['missed']['validated_call'],'1')
            self.assertEqual(rows['missed']['arbitration_status'],'pending')
            reconstruction=targeted_local_reconstruction(bam=root/'evidence/pangenome_reads.bam',region='r3',reference_seq=sequence[2300:3200],outdir=root/'daughter',threads=2)
            self.assertEqual(reconstruction['status'],'reconstructed')
            self.assertGreaterEqual(reconstruction['candidate']['breadth'],.95)
            deleted=sequence[:2300]+sequence[3200:]
            with r1.open('w') as left,r2.open('w') as right:
                for index,start in enumerate(range(0,len(deleted)-400,5)):
                    left.write(f'@deletion{index}\n{deleted[start:start+150]}\n+\n'+('I'*150)+'\n')
                    right.write(f'@deletion{index}\n{deleted[start+250:start+400].translate(complement)[::-1]}\n+\n'+('I'*150)+'\n')
            deletion_bam=root/'deletion.bam'
            map_reads(assembly,str(r1),str(r2),deletion_bam,2,20,root/'deletion.log',retain_ambiguous=True)
            reconstruction=targeted_local_reconstruction(bam=deletion_bam,region='ctg:1801-3700',reference_seq=sequence[2300:3200],outdir=root/'deletion_daughter',threads=2,flank_junction=sequence[1800:2300]+sequence[3200:3700],junction_offset=500)
            self.assertTrue(reconstruction['deletion_spanned'])
