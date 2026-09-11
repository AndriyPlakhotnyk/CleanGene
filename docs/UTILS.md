# CleanGene downstream utilities

`cleangene utils` commands submit SLURM compute jobs by default. The archive commands below also support `--profile local`. The standalone `cleangene-utils` entry point accepts the same subcommands. Results are written under `<run>/results/utils/<analysis-id>` and logs under `<run>/logs/slurm/utils`. Select a run with exactly one of `--run-dir`, `--run ... --analysis-root ...`, or `--latest --analysis-root ...`.

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
  --run <run-id> --analysis-root /path/to/CleanGene \
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

## Wrong-call diagnostics

```bash
cleangene utils diagnose-call \
  --run-dir /path/to/run --organism "Species name" \
  --genes geneA --samples isolate1 isolate2
```

Without `--samples`, CleanGene selects up to `--max-samples` isolates whose initial Panaroo call differs from the read-validated call. The SLURM job identifies the exact BWA-supporting read IDs and sequences, compares raw and processed support, and replays Shovill with retained intermediates. It then checks whether supporting reads survived Shovill's `seqkit` depth sampling and whether they map to the SPAdes graph, SPAdes contigs, replayed Shovill contigs, and original final contigs. Final-contig sequence localization and overlapping Prokka features distinguish assembly, annotation, and Panaroo failures.

Outputs include `diagnostic_summary.tsv`, `supporting_reads.tsv`, `source_artifacts.tsv`, tool versions, mapping artifacts, and per-isolate Shovill replay directories. `--skip-assembly-replay` performs the lighter raw/processed/BWA/final-assembly audit. KMC is reported separately because it estimates k-mer/genome statistics; the direct Shovill read-removal step is depth sampling by `seqkit`.

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

Matrix utilities use `UTILS_CPUS`, `UTILS_MEM`, and `UTILS_TIME`. Read-backed variant analyses use `UTILS_VARIANT_*`; wrong-call replays use `UTILS_DIAGNOSTIC_*`. Account and partition behavior comes from the original run configuration.

If `--cleanup_trimmed_fastq` is passed to `run`/`resume`, or
`CLEANUP_TRIMMED_FASTQ=true`, the final workflow summary replaces retained
fastp outputs with symlinks to the original manifest FASTQs. Read-backed
utilities keep using the paths recorded in isolate QC, but after cleanup they
analyze the original untrimmed reads. The same cleanup can be applied to a
completed older run with `cleangene cleanup --run-dir /path/to/run`.

Assembly compression is controlled during preprocessing with
`COMPRESS_ASSEMBLY_OUTPUTS=off|intermediates|all`, or the matching
`--compress_assembly_outputs` CLI option on `run`/`resume`. The
`intermediates` mode compresses Shovill/SPAdes leftovers such as
`spades.fasta` and `spades.gfa` and does not affect downstream utilities. The
`all` mode also compresses the final `contigs.fa`; CleanGene records
`contigs.fa.gz` in QC and its variant/diagnostic utilities can read that
compressed assembly path. The final summary repeats this operation as a safe,
idempotent sweep so that preprocess tasks completed before a resumed controller
are covered. A report is written to `results/cohort/storage_cleanup.tsv`.

Annotation compression is controlled with
`COMPRESS_ANNOTATION_OUTPUTS=off|nonessential`, or
`--compress_annotation_outputs nonessential`. CleanGene keeps Prokka `.gff`
files uncompressed because Panaroo consumes them, preprocess resume checks look
for them, and variant/diagnostic utilities read neighboring CDS features from
them. The nonessential mode gzips other Prokka outputs such as `.sqn`, `.gbk`,
`.err`, `.ffn`, and `.fna`; the final summary sweep also covers annotations
created before a resumed controller starts.

## Archived alignments and read-derived MSA

```bash
cleangene-utils restore-bam --run-dir /path/to/run --organism "Species name" --samples isolate1 --profile slurm
cleangene-utils inspect-reads --run-dir /path/to/run --organism "Species name" --samples isolate1 --region contig1:100-900 --profile slurm
cleangene-utils evidence-msa --run-dir /path/to/run --organism "Species name" --genes geneA geneB --profile slurm
```

Use `--profile local` for local execution. Omit `--samples` to select all retained
isolates in the organism. `--region` uses each selected isolate's assembly contig
names; omit it to export all alignment records. These operations require completed
CRAM archives from the revised validation workflow. Updating code alone does not
create archives for a completed older run; resume that run to apply current
validation rules and generate archives.

Each operation checks SHA-256 checksums for the CRAM, CRAI, and retained assembly
reference, and writes `alignment_evidence.tsv` with their paths and the reference
checksum. `restore-bam` writes an indexed BAM under each isolate's output directory,
verifies its records against the original archive, and preserves the CRAM.
`inspect-reads` exports `reads.sam` directly from CRAM, optionally by region.

`evidence-msa` restores temporary BAMs, reconstructs consensus with the run's depth,
mapping-quality and base-quality thresholds, lifts known CDS coordinates through
indels, and aligns available gene sequences with MAFFT. Validated discovered CDS
are loaded from their discovery sequences. It writes `<gene>.fasta`,
`<gene>.aligned.fasta`, and `msa_summary.tsv`; unsupported or missing-coordinate
sequences may be absent, so inspect sequence counts. Temporary BAMs are removed
after successful consensus extraction. Retained assembly references and CRAM
files remain available for other read viewers.
