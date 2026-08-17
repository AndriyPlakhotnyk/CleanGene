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

- If `group_id` is omitted, CleanGene uses `organism` as the group label, falling back to `default`.
- Input reads may be supplied as paired FASTQs (`R1`/`R2`) or as an unmapped/raw sequencing BAM (`raw_bam`, `BAM`, or `bam`). Raw BAM reads are extracted with `samtools fastq` before assembly, taxonomy, and BWA validation.
- `READ_TRIMMING_MODE` controls adapter handling: `auto` runs `fastp --detect_adapter_for_pe` when available, `always` requires `fastp`, and `off` leaves reads unchanged.
- A manifest may provide `pangenome_dir`, `panaroo_dir`, or `pangenome` for a group that already has Panaroo output. CleanGene reuses `gene_presence_absence.csv` and the Panaroo sequence files, then runs BWA/SAMtools/BCFtools/minimap2 read validation.
- Kraken2 is an isolate-level taxonomy/contamination gate and requires a site database unless taxonomy is disabled.
- Population-structure inference is not part of CleanGene; its responsibility is the evidence-backed gene matrix and QC.
