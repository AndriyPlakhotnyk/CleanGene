# CleanGene SLURM architecture

The launcher submits one lightweight controller job. The controller polls all running and pending jobs for the user, counts array elements individually, and submits only work that fits below `SLURM_USER_JOB_LIMIT - SLURM_JOB_HEADROOM`.

Arrays use a rolling window rather than a full-chunk barrier. `SLURM_PREPROCESS_MAX_INFLIGHT` and `SLURM_VALIDATION_MAX_INFLIGHT` cap stage occupancy, `SLURM_MAX_PARALLEL` caps execution inside each array, and `SLURM_MAX_OUTSTANDING_CHUNKS` limits active array parents. As array elements complete, the controller refills free capacity without waiting for every element in the older array. QOS submission errors are retried after a fresh all-user occupancy check.

For known organism assignments, the controller orders preprocessing by increasing group size. A group becomes eligible for Panaroo only when every isolate in that group has a successful preprocess marker. Its validation, reduction, and plotting then advance from completion markers while unrelated groups continue preprocessing. Kraken-inferred groups retain the global resolve barrier because their memberships do not exist until Kraken preprocessing completes.

ARC preprocessing uses node-local scratch for generated reads, Kraken reports, Shovill work, assembly, annotation, and logs, then copies final artifacts into the run. Kraken2 uses all allocated CPUs, suppresses per-read classifications unless `KRAKEN2_KEEP_CLASSIFICATIONS=true`, and in `auto` access mode copies the database once per compute node under a lock. If node-local capacity is insufficient, it falls back to memory mapping the shared database.

Worker outputs are checkpointed with JSON completion records. Resubmission is safe after technical failures because workers verify their required final outputs before skipping.
