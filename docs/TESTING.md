# Testing CleanGene

Run the regression suite inside the CleanGene environment:

```bash
bash tests/run_tests.sh
```

External-tool integration tests run when their tools are available. The real
CheckM2 prediction test additionally requires its companion environment and a
completed managed database. A skipped test does not verify that component.

## Installation and updates

The installer supports normal named environments and isolated prefixes. For a
local deployment, run the following twice: the first invocation creates missing
environments; the second exercises updates while preserving local configuration.

```bash
bash scripts/install_or_update.sh --env-root /path/to/test-environments --profile local
bash scripts/install_or_update.sh --env-root /path/to/test-environments --profile local
conda activate /path/to/test-environments/cleangene
```

Omit `--env-root` to use the standard `cleangene` and `cleangene-checkm2` Conda
environments. On ARC use `--profile slurm`; the installer then checks for `sbatch`.

## Real Shovill and CheckM2 pipeline test

With the environment active, use a fresh work directory:

```bash
python scripts/test_shovill_checkm2_e2e.py \
  --work-dir /path/to/e2e-test \
  --checkm2-executable "$CONDA_PREFIX/../cleangene-checkm2/bin/checkm2" \
  --database-root /path/to/shared/checkm2 \
  --threads 2
```

This test uses a 200 kb excerpt of a bundled CheckM2 test genome to generate two
200× paired-read samples with fixed seeds. Shovill exercises genome-size
estimation, subsampling, correction, assembly, and polishing. CleanGene runs real
CheckM2 prediction, annotation, Panaroo, read validation, CRAM archival, and report
generation. It then resumes locally and checks that completed stage markers and
validated calls are reused. `e2e_report.json` is written only after all assertions
pass. Logs and run outputs remain in the work directory for inspection.

The cropped-genome fixture has its completeness exclusion thresholds set to zero;
measured CheckM2 scores remain in the QC report. This is an execution test, not a
benchmark of completeness estimation or gene-call accuracy. An optional `--genome`
selects another input FASTA for the cropped fixture.

The first run may download the full compatible CheckM2 database. Downloads retain
partial bytes for retry, verify the archive checksum, then verify the extracted
DIAMOND database against the SHA-256 recorded in the installed CheckM2 release.
No partially extracted database is exposed as ready. Subsequent runs reuse the
database and runtime verification when their recorded signatures still match.

## Verified local run — 2026-09-11

The pinned environments completed real installation and update checks. The
CheckM2 v3 archive was downloaded from Zenodo with a resumed transfer, and both
archive MD5 and database SHA-256 passed. CheckM2's three bundled test genomes and
its production prediction command completed successfully using companion DIAMOND
2.1.11. All 229 regression tests passed with no skips.

The two-isolate test completed Shovill 1.4.2, CheckM2, Prokka, Panaroo, read
validation, CRAM archival, and summaries: 186 gene clusters, two evaluated
isolates, and two verified archives. Local resume preserved the completed stage
markers and validated matrix. The cropped-fixture limitation above applies.
