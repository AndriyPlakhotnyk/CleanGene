# CleanGene downstream utilities

Every `cleangene utils` command performs lightweight validation on the login node and submits one SLURM compute job. Results are written under `<run>/results/utils/<analysis-id>` and logs under `<run>/logs/slurm/utils`. Select a run with exactly one of `--run-dir`, `--run ... --analysis-root ...`, or `--latest --analysis-root ...`.

## Samples and gene calls

```bash
cleangene utils get-samples \
  --run-dir /path/to/run \
  --organism "Streptococcus gallolyticus" \
  --genes rpoE tetM \
  --match all --status both
```

`gene_presence_absence.tsv` contains the requested genes by isolate. `samples_present.tsv` and `samples_absent.tsv` use either the all-gene or any-gene rule selected with `--match`.

## Differential genes

A selected cohort can be compared with the remaining isolates by using `--samples`, or two explicit cohorts can be supplied with `--cohort-a` and `--cohort-b`. A manifest comparison uses `isolate_id`, `organism`, and `group_id`; `--group-column` changes the grouping column. Exactly two manifest groups are used unless `--group-a-label` and `--group-b-label` select them.

```bash
cleangene utils get-differential-genes \
  --run 260817_120729_cleangene --analysis-root /work/project/CleanGene \
  --manifest cohort.tsv --max-q-value 0.05
```

The comparison reports group counts, prevalence difference, Fisher exact odds ratio and p-value, and Benjamini-Hochberg q-value in `differential_genes.tsv`, with SVG/PNG heatmaps. An organism without cohorts produces a ranked within-organism variable-gene table and heatmap.

## Operons

```bash
cleangene utils get-operon \
  --run-dir /path/to/run --organism "Species name" \
  --operon-name capsule --genes geneA geneB geneC
```

Each distinct binary gene pattern gets an operon ID. Outputs are `operon_calls.tsv`, `operon_summary.tsv`, and SVG/PNG heatmaps. An operon manifest may contain `operon_name`, `genes` (comma/semicolon separated) or `gene`, `organism`, and `isolate_id`. `organism` may be omitted when the supplied sample IDs resolve to one run organism; sample IDs may be omitted to use all isolates in the organism.

## Gene variants

```bash
cleangene utils get-variants \
  --run-dir /path/to/run --organism "Species name" \
  --genes rpoE --min-similarity 95 \
  --flanking-genes 5 --analyze-flanks
```

`--operon /path/to/get-operon-result` selects the genes and cohort from a prior operon analysis. CleanGene recovers Panaroo representative sequences, maps retained paired reads with BWA, calls a haploid BCFtools consensus, and reports identity, coverage, shared bases, SNPs, MNPs, insertions, deletions, truncation bases, assembly coordinates, and neighboring Prokka features. Up to ten neighboring genes per side can be reported. With `--analyze-flanks`, the first observed sequence for each annotated gene and relative flank position becomes a read-validation reference and receives the same variant metrics across the cohort.

Unique qualifying sequences receive stable analysis-local variant IDs. `unique_variants.fasta`, a MAFFT alignment, `unique_variant_alignment.pdf`, `gene_variants.tsv`, and `variant_summary.tsv` link the sequence and tabular results.

## iTOL datasets

```bash
cleangene utils itol \
  --run-dir /path/to/run --organism "Species name" \
  --genes rpoE tetM --color-scheme muted \
  --operon /path/to/operon-result \
  --variants /path/to/variant-result
```

Gene calls are emitted as an iTOL `DATASET_BINARY`; operon and variant IDs are emitted as `DATASET_COLORSTRIP` files. Schemes are `classic`, `muted`, or `custom`. A custom TSV/CSV contains `label` and `color` columns.

## Resources

Matrix utilities use `UTILS_CPUS`, `UTILS_MEM`, and `UTILS_TIME`. Read-backed variant analyses use `UTILS_VARIANT_CPUS`, `UTILS_VARIANT_MEM`, and `UTILS_VARIANT_TIME`. Account and partition behavior comes from the original run configuration.
