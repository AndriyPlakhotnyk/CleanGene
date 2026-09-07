# Locus-aware validation audit

Audited on 2026-09-07 on `feature/locus-aware-gene-validation`.
Fetched origin before editing. The branch already contained the requested
implementation (`8bb67eb`) directly on the latest fetched main (`c21ec0d`), so
this audit continued on that branch. It did not reset the existing work.

The original implementation was partial: secondary-only mappings could disappear
from coverage, 90–95% breadth at high identity was incorrectly called divergent,
consensus indels shifted downstream CDS coordinates, and daughter reconstruction
lacked output directories and name-collated read pairs. These defects are repaired.

## Files changed in this audit

- `src/cleangene/evidence.py`: explicit alignment filtering, batched CDS evidence,
  consensus coordinate liftover and gapped sequence comparison, mapping reuse,
  paralog locus selection, reconstruction repair and deletion junction checks.
- `src/cleangene/workers.py`: arbitration routing and unresolved family handling,
  evidence-version invalidation, reconstructed metric updates, initial-call export.
- `src/cleangene/defaults.py`, `config/cleangene.example.env`: reconstruction
  read/memory bounds and configurable deletion thresholds.
- `tests/test_core.py`, `tests/test_locus_validation.py`: updated coverage contract,
  classification, consensus indels, masks, overlapping CDSs, reconstruction,
  local dispatch, resume, and real synthetic-read integration.
- `README.md`, this report: behavior, configuration, testing, and limitations.

## Biological decisions

Own-assembly CDS evidence takes priority for initial positives. Every assigned
CDS is measured; when a cluster has several copies, the best-supported locus
supplies its aggregate metrics. All assignments remain in `cluster_isolate_loci.tsv`.
Negative calls use competitive pangenome recovery. High-MAPQ primary evidence is
measured separately from family evidence including secondary alignments.
Broad family support with insufficient unique breadth remains `ambiguous_multimap`.

Defaults retain confirmed presence at 95% breadth/identity with mean depth >=5;
divergent variants require >=90% breadth and identity from 90% to below 95%.
High-identity breadth from 70% to below 95% is `possible_truncation`; initial
positives remain provisionally present. `partial_homolog` retains an intact-gene
call of zero and is never presented as a demonstrated deletion.

Initial positive `not_detected` cases enter arbitration and retain the initial
binary call if unresolved. Strong negative-to-positive recoveries and ambiguous
cases also enter the shared local/Slurm daughter stage. Family-only arbitration
keeps `validated_call` empty; the compatibility binary matrix carries forward the
initial call. An empty validated call is not biological absence.

Daughter reconstruction retains mates, collates names before FASTQ export,
combines singleton sources, and runs bounded SPAdes. A deletion requires a locally
reconstructed sequence with a contiguous alignment across both flanks. An insertion
or deletion separating the two alignment anchors does not prove the junction.
An ambiguous-family reconstruction alone cannot resolve an exact cluster.

The command contracts were checked against the official
[SAMtools coverage documentation](https://www.htslib.org/doc/1.20/samtools-coverage.html),
[SAMtools depth documentation](https://www.htslib.org/doc/1.20/samtools-depth.html),
and [BCFtools consensus documentation](https://samtools.github.io/bcftools/bcftools).

## Outputs and configuration

Existing binary filenames remain unchanged. `read_validation_metrics.tsv` now
includes the previously omitted `initial_call`. Rich outputs include evidence and
sequence-resolution states, validated calls, call source, breadth, depth,
normalized depth, reconstructed identity/length, reference length, ORF status,
unique and ambiguous alignment counts, scaffold/CDS coordinates, contig-edge
status, and arbitration status/reason.

This audit adds `evidence_version`, `family_breadth`, and
`reconstructed_coverage`. The latter is observed consensus bases/reference length,
capped at one; it is distinct from read breadth. Alignment counts are not counts
of distinct DNA molecules. For a locus-positive row, `family_breadth` currently
reports the high-MAPQ locus breadth; ambiguous alignment counts are retained
separately. For recovery rows, family breadth includes ambiguous mappings.

New configuration variables:

| Variable | Default | Meaning |
| --- | --- | --- |
| `READ_VALIDATION_ARBITRATION_MAX_READS` | 100000 | Maximum recruited read names per case; larger cases stay unresolved |
| `READ_VALIDATION_ARBITRATION_MEMORY_GB` | 12 | SPAdes memory limit, below the default 16G job request |
| `READ_VALIDATION_DELETION_MIN_IDENTITY` | 0.95 | Minimum reconstructed junction identity |
| `READ_VALIDATION_DELETION_MIN_ANCHOR` | 50 | Contiguous aligned bases required on each side |

Existing controls include the six breadth/depth/identity thresholds,
`READ_VALIDATION_MIN_MAPQ`, `BASEQUAL`, `READ_VALIDATION_FLANK_LENGTH`,
`READ_VALIDATION_ARBITRATION_MAX_CASES`, `SLURM_ARBITRATION_MAX_INFLIGHT`,
and `ARBITRATION_CPUS`, `ARBITRATION_MEM`, `ARBITRATION_TIME`.

## Computational impact and resume

Two mappings per isolate are retained: own assembly and pangenome. CDS measurement
uses two streamed BAM passes rather than subprocesses per CDS. Ambiguous
alignments and unmapped mates increase BAM size. Reconstruction is restricted to
discordant cases and bounded by case count, recruited names, and SPAdes memory.
Local and Slurm modes dispatch the same biological workers.

Mapping completion signatures include input paths, sizes, modification times,
MAPQ/retention options, and evidence version. Successful BAM/index pairs with
matching signatures are reused. Resume maintenance invalidates older evidence
and downstream completion markers while preserving unrelated upstream stages.
Consensus files are rebuilt after an interrupted validation; completed isolate
stages retain the existing workflow's marker behavior.

## Validation and remaining limitations

The full suite passes: 197 tests run, one skipped because the real CheckM2
runtime/database test is unavailable. The real synthetic-read integration uses
BWA, SAMtools, BCFtools, minimap2, and SPAdes from the installed CleanGene
environment. It verifies own-locus presence despite competitive ambiguity,
family-only homolog evidence, annotation-missed recovery, daughter assembly,
and a reconstructed deletion junction. Existing tests cover Slurm arrays,
compression, assembly/annotation, QC, resume, and unrelated workflow behavior.

ARC's production manifest, scheduler, and full cohort were not exercised.
The following are explicit limits, not resolved biological conclusions:

- Automatic transfer of conserved flanking loci from other isolates to an
  initial negative is not implemented. Junction testing needs sample coordinates.
- Exact allele discrimination among homologous reconstructions remains unresolved;
  there is no cohort-wide competing-allele arbitration algorithm.
- Own-locus identity is relative to the sample CDS; recovery identity is relative
  to the pangenome reference. These are different biological comparisons.
- Normalized depth uses a length-weighted median of assembly-contig depths as a
  chromosomal proxy, without chromosome/plasmid labeling.
- ORF checks are basic length/stop checks; masked sequence is unresolved.
- Reconstruction uses minimap2 asm5 and can miss short or very divergent genes.
  Global consensus comparisons exceeding two million residual alignment cells
  after common-end trimming remain unresolved to bound computation.
- Preparation still materializes cohort tables and locus mappings; tens-of-thousands
  scaling has not been benchmarked. No throughput claim is made for that size.
- The run's existing completion markers do not generally fingerprint every
  configuration change; use a new run when changing biological thresholds.

## Commands

From this checkout, with the CleanGene environment activated:

```bash
bash tests/run_tests.sh
PYTHONPATH=src python -m unittest discover -s tests -p test_locus_validation.py -v
```

Local integration (replace the manifest path with a validated cohort containing
at least two isolates per organism; use a fresh analysis root):

```bash
PYTHONPATH=src python -m cleangene run \
  --profile local \
  --manifest /path/to/cohort.manifest.tsv \
  --analysis-root /tmp/cleangene-locus-local \
  --config config/cleangene.example.env \
  --assembler spades \
  --compress-assembly-outputs intermediates \
  --compress-annotation-outputs nonessential
```

On ARC, install this branch into the normal CleanGene environment and run the
original regression command with a fresh run ID. This command is documented,
not submitted by this audit:

```bash
CLEANGENE_ROOT="$PWD"  # Run from your CleanGene checkout on ARC.
cleangene run \
  --profile slurm \
  --manifest "${CLEANGENE_ROOT}/input/arc_GDS_test.manifest.tsv" \
  --analysis-root "${CLEANGENE_ROOT}" \
  --config "${CLEANGENE_ROOT}/config/cleangene.arc.no-checkm2.env" \
  --assembler spades \
  --compress-assembly-outputs intermediates \
  --compress-annotation-outputs nonessential
```

Add `--dry-run` to inspect submission construction; remove it to execute. Inspect
`gene_call_evidence.long.tsv`, arbitration cases, and controller logs, alongside
pipeline completion. Binary agreement alone is not a scientific acceptance test.
