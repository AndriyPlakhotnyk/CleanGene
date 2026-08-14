# Validation status

## Automated checks included

- pure-Python unit tests for assembly metrics, Panaroo call normalization, validation-gene selection, and Kraken contamination parsing;
- SLURM dry-run generation exercises the complete dependency chain without submitting jobs;
- Python modules are compile-checked before release packaging.

Run:

```bash
./tests/run_tests.sh
./cleangene run --manifest input/cohort.tsv --config config/cleangene.example.env --analysis-root /path/to/work --profile slurm --dry-run
```

## Required site validation

This source release has not executed Shovill, Prokka, Panaroo, BWA, SAMtools, BCFtools, minimap2, and Kraken2 end-to-end against your production cluster/data from this build environment. Before a 30,000-isolate production run, execute a representative real cohort and compare CleanGene calls against the corresponding ImaGene outputs.

## Deliberate v0.1 boundaries

- `group_id` is supplied by the user; CleanGene does not infer species/lineage groups.
- Kraken2 is an isolate-level taxonomy/contamination gate and requires a site database unless taxonomy is disabled.
- Population-structure inference is not part of CleanGene; its responsibility is the evidence-backed gene matrix and QC.
