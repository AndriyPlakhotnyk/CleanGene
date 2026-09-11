import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cleangene.defaults import DEFAULTS
from cleangene.evidence import classify_gene_evidence
from cleangene.util import atomic_json, read_tsv, write_tsv
from cleangene.validation_summary import (SUMMARY_FILES, PLOT_FILES, decision_metric_membership,
    write_validation_summary, combine_validation_summaries)
from cleangene.workers import invalidate_missing_validation_reports


class ValidationSummaryTests(unittest.TestCase):
    def fixture(self, root, old, new, isolates=('a','b'), records=None):
        initial=root/'initial.tsv'; final=root/'validated_gene_presence_absence.binary.tsv'; evidence=root/'gene_call_evidence.long.tsv'
        write_tsv(initial,['Gene',*isolates],old); write_tsv(final,['Gene',*isolates],new)
        if records is None:
            records=[]
            for before,after in zip(old,new):
                for index,iso in enumerate(isolates,1):
                    state='confirmed_present' if after[index] else 'partial_homolog'
                    records.append({'Gene':before[0],'isolate_id':iso,'initial_call':before[index],'final_call':after[index],
                                    'evidence_state':state,'decision_reason':state.replace('_',' '),'decision_metrics':'breadth;identity;mean_depth'})
        write_tsv(evidence,['Gene','isolate_id','initial_call','final_call','evidence_state','decision_reason','decision_metrics'],records)
        return initial,final,evidence

    def test_both_directions_kept_zeros_and_all_gene_denominator(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            paths=self.fixture(root,[['g1',0,1],['g2',1,0],['g3',0,0],['g4',1,1]],
                               [['g1',1,1],['g2',0,0],['g3',0,0],['g4',1,0]])
            totals=write_validation_summary(*paths,root,'group')
            self.assertEqual([totals[k] for k in ('total_calls','kept','changed','changed_0_to_1','changed_1_to_0')],[8,5,3,1,2])
            self.assertEqual(totals['changed_pct'],37.5)
            per=read_tsv(root/'gene_call_summary.per_isolate.tsv')
            self.assertEqual([r['changed_pct'] for r in per],['50.0','25.0'])
            self.assertEqual([r['kept'] for r in per],['2','3'])
            self.assertIn('Validation: PASS',(root/'summary_statistics.txt').read_text())
            self.assertEqual(sum(int(row['n_calls']) for row in read_tsv(root/'decision_reason_upset.tsv')),3)

    def test_prevalence_boundaries_and_absent_clusters_from_final_matrix(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); isolates=tuple(f'i{x}' for x in range(100))
            final=[[f'g{n}',*([1]*n+[0]*(100-n))] for n in (100,99,95,94,15,14,0)]
            initial=[[r[0],*([1]*100)] for r in final]
            paths=self.fixture(root,initial,final,isolates)
            write_validation_summary(*paths,root,'g'); report=(root/'summary_statistics.txt').read_text()
            for line in ('Core genes\t(99% <= isolates <= 100%)\t2','Soft core genes\t(95% <= isolates < 99%)\t1',
                         'Shell genes\t(15% <= isolates < 95%)\t2','Cloud genes\t(0% <= isolates < 15%)\t2',
                         'Absent in all isolates (included in cloud)\t1'):
                self.assertIn(line,report)

    def test_evidence_mismatch_fails_without_success_report(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); paths=self.fixture(root,[['g',0,0]],[['g',1,0]])
            text=paths[2].read_text().replace('g\ta\t0\t1','g\ta\t0\t0');paths[2].write_text(text)
            (root/'summary_statistics.txt').write_text('old success')
            with self.assertRaisesRegex(ValueError,'contradicts'): write_validation_summary(*paths,root,'g')
            self.assertFalse((root/'summary_statistics.txt').exists())

    def test_missing_and_duplicate_matrix_rows_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); paths=self.fixture(root,[['g',0,1]],[['g',0,1],['g',0,1]])
            with self.assertRaisesRegex(ValueError,'Matrix gene rows'): write_validation_summary(*paths,root,'g')

    def test_empty_matrix_and_no_changes_have_valid_reports_and_plot(self):
        from cleangene.plotting import plot_decision_upset
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);paths=self.fixture(root,[],[])
            totals=write_validation_summary(*paths,root,'empty')
            self.assertEqual(totals['total_calls'],0); self.assertEqual(totals['kept_pct'],0)
            plot_decision_upset(root/'decision_reason_upset.tsv',root,'empty')
            self.assertTrue(all((root/name).stat().st_size>0 for name in PLOT_FILES))

    def test_legacy_membership_is_labeled_not_fabricated(self):
        metrics,source=decision_metric_membership({'evidence_state':'confirmed_present'})
        self.assertEqual(source,'inferred_from_state'); self.assertIn('identity',metrics)
        metrics,source=decision_metric_membership({'decision_reason':'unknown old decision'})
        self.assertEqual(source,'unavailable'); self.assertEqual(metrics,('unrecorded_evidence',))
        result=classify_gene_evidence(mapped_reads=0,breadth=0,mean_depth=0,identity=None)
        self.assertEqual(result['decision_metrics'],'mapped_reads;breadth')
        self.assertEqual(decision_metric_membership(result)[1],'recorded')

    def test_cohort_percentages_are_weighted_by_calls(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); dirs=[root/'g1',root/'g2']
            for out in dirs: out.mkdir()
            write_validation_summary(*self.fixture(dirs[0],[['g',0,0]],[['g',1,1]]),dirs[0],'g1')
            old=[['g1',1,1],['g2',1,1],['g3',0,0]]
            write_validation_summary(*self.fixture(dirs[1],old,old),dirs[1],'g2')
            combine_validation_summaries(dirs,root/'cohort')
            row=read_tsv(root/'cohort/gene_call_summary.tsv')[0]
            self.assertEqual(row['total_calls'],'8');self.assertEqual(row['changed_pct'],'25.0')
            self.assertEqual(len(read_tsv(root/'cohort/gene_call_summary.per_isolate.tsv')),4)

    def test_resume_backfills_only_reporting_stages(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); out=root/'results/groups/g/03_read_validation'; out.mkdir(parents=True)
            (out/'validated_gene_presence_absence.binary.tsv').write_text('Gene\ti\ng\t1\n')
            for stage in ('validate','arbitrate','reduce','plot'): atomic_json(root/f'state/{stage}/g.done.json',{})
            atomic_json(root/'state/summary.done.json',{})
            self.assertEqual(invalidate_missing_validation_reports(root),1)
            for stage in ('validate','arbitrate'): self.assertTrue((root/f'state/{stage}/g.done.json').exists())
            for stage in ('reduce','plot'): self.assertFalse((root/f'state/{stage}/g.done.json').exists())
            for name in (*SUMMARY_FILES,*PLOT_FILES): (out/name).touch()
            self.assertEqual(invalidate_missing_validation_reports(root),0)

    def test_compression_defaults_and_explicit_off(self):
        from cleangene.config import read_env
        self.assertEqual(DEFAULTS['COMPRESS_ASSEMBLY_OUTPUTS'],'intermediates')
        self.assertEqual(DEFAULTS['COMPRESS_ANNOTATION_OUTPUTS'],'nonessential')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'config.env';p.write_text('COMPRESS_ASSEMBLY_OUTPUTS=off\nCOMPRESS_ANNOTATION_OUTPUTS=off\n')
            cfg=read_env(p)
            self.assertEqual(cfg['COMPRESS_ASSEMBLY_OUTPUTS'],'off')
            self.assertEqual(cfg['COMPRESS_ANNOTATION_OUTPUTS'],'off')

    def test_plot_contains_directions_reasons_and_metric_matrix(self):
        from cleangene.plotting import plot_decision_upset
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);paths=self.fixture(root,[['g',0,1]],[['g',1,0]])
            write_validation_summary(*paths,root,'example');plot_decision_upset(root/'decision_reason_upset.tsv',root,'example',max_columns=1)
            svg=(root/'decision_reason_upset.svg').read_text()
            for label in ('confirmed present','identity','Call change','decision_reason'):
                self.assertIn(label,svg)
            self.assertTrue((root/'decision_reason_upset_page_02.svg').is_file())

    def test_reduce_and_plot_publish_reports_and_backfill_without_read_work(self):
        from cleangene.workers import reduce_group, plot_group
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); group='g'; root=run/'results/groups/g'; out=root/'03_read_validation'
            atomic_json(run/'provenance/resolved_config.json',{})
            write_tsv(root/'02_pangenome/initial_calls/gene_presence_absence.binary.tsv', ['Gene','a','b'],
                      [['g1',0,1],['g2',1,0]])
            fields=['Gene','validated_call','evidence_state','decision_reason','decision_metrics','final_call_source']
            write_tsv(out/'evidence/a/metrics.tsv',fields,[
                ['g1',1,'confirmed_present','confirmed present','breadth;identity;mean_depth','pangenome_read_recovery'],
                ['g2',0,'partial_homolog','partial homolog','breadth;identity;mean_depth','targeted_local_reconstruction']])
            write_tsv(out/'evidence/b/metrics.tsv',fields,[
                ['g1','','ambiguous_multimap','ambiguous multimap','unique_mapped_reads;ambiguous_mapped_reads','initial_call_unresolved']])
            retained=[{'isolate_id':iso} for iso in ('a','b')]
            with patch('cleangene.workers.task_row',return_value={'group_id':group}), patch('cleangene.workers.retained_rows',return_value=retained), patch('cleangene.workers.validate_isolate',side_effect=AssertionError('unnecessary read validation')):
                reduce_group(run,0); plot_group(run,0)
                self.assertTrue(all((out/name).is_file() for name in (*SUMMARY_FILES,*PLOT_FILES)))
                summary=read_tsv(out/'gene_call_summary.tsv')[0]
                self.assertEqual(summary['changed'],'2');self.assertEqual(summary['kept'],'2')
                self.assertEqual(read_tsv(root/'cleaned_pangenome.tsv')[0]['b'],'1')
                modified=(out/'summary_statistics.txt').stat().st_mtime_ns
                reduce_group(run,0)
                self.assertEqual(modified,(out/'summary_statistics.txt').stat().st_mtime_ns)
                (out/'summary_statistics.txt').unlink()
                self.assertEqual(invalidate_missing_validation_reports(run),1)
                reduce_group(run,0); plot_group(run,0)
                self.assertEqual(read_tsv(out/'gene_call_summary.tsv')[0],summary)
