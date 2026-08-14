# CleanGene SLURM architecture

The launcher submits coarse arrays, not one job per gene. Isolate preprocessing and read validation are arrays; Panaroo, validation preparation and reduction are arrays over groups. Each stage uses `afterok` dependencies.

`SLURM_MAX_PARALLEL` caps array concurrency. Panaroo has separate CPU/memory/time controls because its memory behavior differs substantially from per-isolate mapping.

Worker outputs are checkpointed with JSON completion records. Resubmission is safe after technical failures because workers verify their required final outputs before skipping.
