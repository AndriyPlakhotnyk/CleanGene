# ARC handoff for local ML agents

This note describes the CleanGene checkout and ARC workflow as of commit
`363ad36` (`Add resistance operon analysis workflow`). The preceding commit
`17805ca` added bounded recovery for Slurm preprocess arrays that reach
`TIMEOUT`.

## Checkout and runtime

- The production checkout is a Git clone on ARC project/shared storage. Keep
  the repository path, Conda environments, databases, run directories and
  input reads on project/shared storage; use node-local scratch only for
  per-task temporary files.
- The tracked ARC template is `config/cleangene.arc.env`. The machine-specific
  `config/cleangene.arc.local.env` is intentionally local and must retain the
  account, partition, database roots and site paths.
- The primary environment is `cleangene`; CheckM2 uses the companion
  `cleangene-checkm2` environment. The shared CheckM2 database is recorded in
  the resolved run configuration and guarded by a shared verification marker
  and lock. Do not delete or replace that database while a run is active.
- A run stores its immutable resolved configuration in
  `provenance/resolved_config.json`, task lists in `state/`, Slurm logs in
  `logs/slurm/`, and developer timing reports in `logs/developer_final_report/`.
  Resume from the run directory so successful markers and outputs are reused.

## ARC scheduler assumptions

The current ARC template uses these defaults; local configuration may override
them for the allocation:

| Setting | Default | Meaning |
| --- | ---: | --- |
| `SLURM_USER_JOB_LIMIT` | 2000 | User-wide running/pending array elements allowed by the controller. |
| `SLURM_JOB_HEADROOM` | 10 | Capacity reserved below that limit. |
| `SLURM_MAX_PARALLEL` | 400 | Concurrent elements within one array parent. |
| `SLURM_PREPROCESS_MAX_INFLIGHT` | 400 | Preprocess occupancy cap. |
| `SLURM_ARRAY_CHUNK_SIZE` | 500 | Maximum task indices in one submitted chunk. |
| `SLURM_MAX_OUTSTANDING_CHUNKS` | 8 | Active array-parent cap per stage. |
| `SLURM_POLL_SECONDS` | 20 | Scheduler polling/refill interval. |
| `SLURM_CONTROLLER_REPORT_INTERVAL_SECONDS` | 120 | Controller progress-report interval. |
| `SLURM_PREPROCESS_TIME` | 24:00:00 | Preprocess wall-time request. |
| `SLURM_TIMEOUT_RETRIES` | 1 | Bounded retry count for missing work after a Slurm `TIMEOUT`. |

Seeing 400/2000 is expected: 400 is the preprocess stage cap, while 2000 is
the user-wide queue limit. The controller refills arrays as elements finish.
Array elements share a parent ID in Slurm; CleanGene tracks task IDs and output
markers separately.

When a parent reaches `TIMEOUT`, the controller reconciles completed outputs,
removes the expired parent, and resubmits only indices without successful
markers. `FAILED`, `CANCELLED`, and out-of-memory states remain fatal. A resume
is therefore safe after a controller failure once no old CleanGene jobs for
that run remain active.

## Latest resistance workflow

`--resistance-operon` is optional and disabled by default. It adds AMRFinderPlus
scan/read-support work after the ordinary pangenome reduction and plots, then
merges resistance loci, similarity bins, exact variants, alignments and figures
before final summary/archiving. It requires a versioned shared `AMRFINDER_DB`
and the updated primary environment. Outputs are under
`results/resistance_analysis/`; definitions and software/database versions are
retained in its provenance directory. The default origin is `unknown` unless
reviewed contig-origin evidence is supplied.

The supported provisioning path is `scripts/bootstrap_arc_latest.sh`, after
loading the ARC Git and Miniforge/Mamba modules. It fast-forwards or clones the
checkout, updates both environments, provisions AMRFinderPlus, writes the
versioned database path into the local config, and runs `cleangene doctor`.

## Commands used on ARC

For a normal existing checkout:

```bash
module load git
module load miniforge3                 # use the module name supplied by ARC
cd <ARC_PROJECT>/pipelines/CleanGene
git status --short                     # preserve local changes before pulling
git pull --ff-only origin main
bash scripts/install_or_update.sh --profile slurm
mamba run -n cleangene cleangene doctor \
  --profile slurm --config config/cleangene.arc.local.env
```

For a clean/shared installation or the resistance workflow, use:

```bash
bash scripts/bootstrap_arc_latest.sh \
  <ARC_PROJECT>/pipelines/CleanGene \
  <ARC_PROJECT>/software/conda \
  <ARC_PROJECT>/databases/amrfinderplus main
```

To resume the interrupted run, substitute its actual run directory:

```bash
mamba run -n cleangene cleangene resume \
  --run-dir <ARC_PROJECT>/pipelines/CleanGene/runs/<RUN_ID> \
  --config config/cleangene.arc.local.env \
  --skip-downsampling
```

Add `--resistance-operon` when adding the resistance stage to an existing run.
Add `--ignore-checkm2` only when the analysis is intentionally configured to
omit CheckM2; it changes the scientific QC record and should be explicit.

## Guidance for future local ML changes

1. Read the resolved run configuration and completion markers before proposing
   work. Do not infer completion from `squeue` alone.
2. Preserve task indices and marker semantics. A worker may be safely rerun only
   after checking its required outputs and input/config provenance.
3. Keep Slurm submission, local execution and resume paths aligned. Test a
   scheduler dry run before submitting a large array.
4. Use `bash tests/run_tests.sh` in the `cleangene` environment. The resistance
   real-tool test is `scripts/test_resistance_e2e.py`; it is a synthetic fixture,
   not an ARC-scale performance benchmark.
5. Treat CheckM2 and AMRFinderPlus databases as shared, versioned resources.
   Diagnose executable/database/runtime mismatches before changing worker logic.
6. Record controller reports, developer timing reports and failed task-log
   excerpts in any handoff. They are the evidence for queue counts, completed
   samples, stage runtimes and resume decisions.
