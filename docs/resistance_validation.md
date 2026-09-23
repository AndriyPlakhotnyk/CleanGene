# Resistance analysis verification — 2026-09-22

- Regression suite: **269 tests passed** (`tests/run_tests.sh`, CleanGene Conda
  environment; local multiprocessing enabled for the existing CheckM2 test).
- Fresh real-tool integration test: **PASS**, including a second pipeline
  invocation that reused both isolate scans and unchanged reports.
- Tested path: generated paired FASTQs → fastp → Shovill → Prokka → Panaroo →
  read validation/arbitration → AMRFinderPlus → locus extraction/read evidence →
  95% bins → MAFFT → figures/tables → final archiving.
- `SKIP_DOWNSAMPLING=true`, `CHECKM2_MODE=off`, `RESISTANCE_OPERON=true` were
  verified in resolved configuration. Shovill logs confirm no depth reduction.
- Two in-silico isolates with a public tet(M) reference yielded two exact sequence
  variants in one similarity bin. The expected SNP and 3-base indel were detected;
  both loci passed the read-support breadth threshold. Origin stayed unknown.
- AMRFinderPlus **4.2.7**, database **2026-08-07.1** were downloaded into `/tmp`
  for testing; the user's installed CleanGene environment was not modified.
- Slurm launcher dry run: **PASS**. It generated an `sbatch` controller command
  and persisted the requested flags. A 999-isolate scheduling test verifies
  per-isolate indices and resistance resource settings. No live ARC submission
  or thousand-isolate performance test was performed.
- Figures were inspected visually: SNP/indel overview, base-level differences,
  and gene order/orientation diagrams. Shell syntax and `git diff --check` pass.

The fresh test report and outputs are available in the working checkout under:

```text
test_data/resistance-e2e-verified/e2e_report.json
test_data/resistance-e2e-verified/runs/e2e/results/resistance_analysis/
```

`test_data/` is ignored by Git. The reproducible test driver is
[`scripts/test_resistance_e2e.py`](../scripts/test_resistance_e2e.py).
The small computational fixture is not a biological E. faecium benchmark and
uses taxonomic classification off. It does not validate Kraken2 classification
accuracy or performance on the forthcoming cohort. See
[methods and submission instructions](resistance_operons.md).
