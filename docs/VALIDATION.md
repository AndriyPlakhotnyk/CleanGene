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

## Read-validation call states

CleanGene writes `results/cohort/validation_decision_logic.tsv` and per-gene evidence files with the same state model:

| State | Criteria | Final call behavior |
| --- | --- | --- |
| `not_detected` | `mapped_reads=0` or `breadth=0` | Set to `0` |
| `low_depth` | some coverage but mean depth below `READ_VALIDATION_MIN_MEAN_DEPTH` | Preserve initial pangenome call |
| `partial_coverage` | depth passes but breadth below `READ_VALIDATION_MIN_BREADTH` | Set to `0` |
| `identity_unresolved` | breadth/depth pass but sequence identity cannot be measured | Preserve initial pangenome call |
| `divergent` | breadth/depth pass but identity below `READ_VALIDATION_MIN_IDENTITY` | Set to `0` |
| `confirmed_present` | breadth, depth, and identity pass | Set to `1` |
| `not_tested_carried_forward` | gene was not selected for read validation | Preserve initial pangenome call |

Resume invalidates legacy metrics where covered genes were serialized as `identity=0` with no aligned positions, then reruns only the affected validation/reduce/plot/summary stages.

## Deliberate v0.1 boundaries

- If `group_id` is omitted, CleanGene uses `organism` as the group label, falling back to `default`.
- Input reads may be supplied as paired FASTQs (`R1`/`R2`) or as one paired unmapped sequencing BAM (`raw_bam`; aliases `ubam`, `uBAM`, `unaligned_bam`, `BAM`, and `bam`). CleanGene requires balanced unmapped READ1/READ2 records, name-collates the uBAM, and extracts complete pairs with `samtools fastq` before assembly, taxonomy, and BWA validation. A row cannot provide both input forms.
- `READ_TRIMMING_MODE` controls adapter handling: `auto` runs `fastp --detect_adapter_for_pe` when available, `always` requires `fastp`, and `off` leaves reads unchanged.
- Skani is not part of the CleanGene pipeline. If `group_id` is absent but `organism` is supplied, `organism` becomes the pangenome group. If both are absent, CleanGene runs Kraken2 and groups isolates by inferred top species.
- With `TAXONOMY_MODE=auto`, Kraken2 runs only when needed for missing organism/group IDs. Use `TAXONOMY_MODE=kraken2` to force taxonomic QC for already grouped manifests.
- A manifest may provide `pangenome_dir`, `panaroo_dir`, or `pangenome` for a group that already has Panaroo output. CleanGene reuses `gene_presence_absence.csv` and the Panaroo sequence files, then runs BWA/SAMtools/BCFtools/minimap2 read validation.
- Kraken2 is an isolate-level taxonomy/contamination gate and requires a site database unless taxonomy is disabled.
- Population-structure inference is not part of CleanGene; its responsibility is the evidence-backed gene matrix and QC.
