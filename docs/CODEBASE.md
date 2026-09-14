# Codebase and verified execution baseline

CleanGene works end-to-end locally. The retained real-tool integration reports
verify the complete pipeline and completed-stage reuse, not just mocked workers.
The September 2026 baseline includes:

| Fixture | Result | Gene clusters | CheckM2 isolates | Verified CRAM archives | Resume reuse |
| --- | --- | --- | --- | --- | --- |
| Normal Shovill processing | PASS | 186 | 2 | 2 | Yes |
| Shovill with downsampling disabled | PASS | 185 | 2 | 4 | Yes |

The second fixture also restored four BAMs from their archives. Both use a 200 kb
cropped genome with completeness exclusion thresholds set to zero. These results
verify execution, not biological accuracy or cohort-scale performance. The
original reports and sequencing artifacts are local, ignored data; reproduce them
using [TESTING.md](TESTING.md). Do not substitute unit-test success for a fresh
real-tool end-to-end result.

## Organization

- `src/cleangene/cli.py`: run creation, configuration, local stage dispatch and
  Slurm submission; `utils_cli.py` exposes downstream utilities.
- `workers.py`: shared preprocessing, grouping, Panaroo, validation, arbitration,
  reduction, plotting and summary workers; also the rolling Slurm controller.
- `slurm.py`, `task_store.py`, `completion.py`: scheduling, stable task lookup,
  output validation and resume reconciliation.
- `evidence.py`, `discovery.py`, `pangenome.py`: read evidence classification,
  reconstruction/discovery and matrix/reference handling.
- `qc.py`, `checkm2.py`, `kraken.py`: sample QC and external tool integration.
- `plotting.py`, `validation_summary.py`: matrix figures and reconciled reports.
- `alignment_archive.py`, `final_archives.py`, `archive_utils.py`: verified CRAM
  archival and restoration; `downstream.py` implements downstream analyses.
- `config/` contains deployment templates; `scripts/` contains installation and
  real-tool integration helpers; `tests/` contains regression/integration checks.

## Transfer to Slurm

Local execution and Slurm dispatch the same stage workers and biological
configuration. Local execution is sequential; Slurm uses dependencies, arrays,
resource tiers and a rolling controller. Before deployment:

1. Install the same committed checkout with
   `bash scripts/install_or_update.sh --profile slurm` and activate `cleangene`.
2. Set account, partition and stage resources in a private copy of
   `config/cleangene.arc.env`. Put input reads, databases and the analysis root on
   storage visible to all workers. Paths and symlink targets from a workstation
   are not automatically portable.
3. Run `cleangene doctor --profile slurm --config CONFIG`, then submit a small
   manifest with `cleangene run --profile slurm --manifest MANIFEST
   --analysis-root ROOT --config CONFIG`. Preserve biological settings, including
   `--skip-downsampling` when reproducing that baseline.
4. Verify final matrices, evidence, before/after figures and summary markers;
   exercise `cleangene resume --run-dir RUN --config CONFIG` and confirm completed
   biological stages remain reused.

Slurm orchestration has regression coverage; an actual cluster execution is a
separate deployment check and has not been newly verified on this workstation.
Prefer a new cluster run with relocated manifest paths to copying a completed
local run containing absolute paths. Never update worker code while jobs using
that checkout are active. Keep private manifests, databases, generated results,
backup source files and machine-specific AI notes out of commits.
