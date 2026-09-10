# CleanGene

**From bacterial sequencing reads to pangenomes you can inspect.**

CleanGene combines read and assembly quality control, annotation, Panaroo pangenome
construction, and read-backed gene validation. It produces a binary gene
presence/absence matrix alongside evidence that distinguishes supported genes,
partial homologs, divergent variants, ambiguous mappings, and unresolved calls.

Run cohorts through Slurm or use local execution for smaller analyses. Both modes
produce the same types of gene evidence.

| Start with | Analyze | Receive |
| --- | --- | --- |
| Paired FASTQs or unmapped BAMs | Read QC, assembly, annotation, and Panaroo | A cleaned gene presence/absence matrix |
| Existing assemblies or Panaroo results | Locus support and targeted reconstruction | Evidence tables explaining each gene call |

[Pipeline](#pipeline) · [Gene decisions](#gene-presence-and-absence-decisions) ·
[Installation](#installation) · [Quick start](#quick-start) · [Arguments](#command-line-arguments) ·
[Outputs](#outputs)

## Pipeline

```mermaid
flowchart TD
    A["Manifest: paired FASTQs or uBAM"] --> B{"Execution profile"}
    B -->|Slurm| C["Submit controller and stage arrays"]
    B -->|Local| D["Run stages locally"]
    C --> E["Preflight and required database setup"]
    D --> E
    E --> F["Read QC and processing; Kraken2 taxonomy"]
    F --> G["Assembly and annotation; CheckM2 when enabled"]
    G --> H{"Isolate QC"}
    H -->|FAIL| X["Record exclusion and QC evidence"]
    H -->|PASS or WARNING| I["Resolve retained analysis groups"]
    I --> J["Run Panaroo or load supplied Panaroo output"]
    J --> K["Prepare gene references and sample CDS coordinates"]
    K --> L["Map reads to own assembly and pangenome references"]
    L --> M["Classify gene evidence"]
    M --> N{"Discordance or unresolved assignment?"}
    N -->|Yes| O["Bounded targeted reconstruction and arbitration"]
    N -->|No| P["Combine gene calls"]
    O --> P
    P --> Q["Binary matrices, evidence tables, and cohort summaries"]
```

The diagram shows stage dependencies. Slurm may overlap independent samples and
groups. Existing artifacts can bypass their corresponding preparation stages.

Isolates receive `PASS`, `WARNING`, or `FAIL` QC status. Warnings remain eligible
for downstream analysis; failures are excluded. QC includes taxonomy, read
quality, sequencing coverage, assembly quality, annotation success, and CheckM2
when enabled. Thresholds can be configured globally or through QC profiles and
per-isolate manifest overrides.

## Gene presence and absence decisions

Panaroo supplies the initial calls. Initial positives are evaluated at their
sample-specific CDS coordinates using reads mapped to the isolate's own assembly.
Initial negatives are searched against pangenome references to recover genes
missed by assembly or annotation. Missing CDS coordinates also require the
reference-search fallback.

```mermaid
flowchart TD
    A{"Initial Panaroo call"} -->|Present| B["Own-assembly CDS support"]
    A -->|Absent| C["Pangenome reference search"]
    B --> D["Evaluate breadth, depth, consensus identity, and mapping ambiguity"]
    C --> D
    D --> E{"Family supported but exact assignment unresolved?"}
    E -->|Yes| F["ambiguous_multimap"]
    E -->|No| G["Classify sequence evidence using configured thresholds"]
    G --> H["confirmed_present or divergent_variant: call 1"]
    G --> I["possible_truncation: retain initial positive provisionally"]
    G --> J["partial_homolog: intact-gene call 0"]
    G --> K["not_detected or insufficient_evidence"]
    F --> L["Arbitrate discordant and unresolved cases"]
    I --> L
    J -->|Initial positive| L
    K -->|Initial positive| L
    H -->|Recovered initial negative| L
    L --> M{"Targeted reconstruction outcome"}
    M -->|Supported deletion junction| N["confirmed_absent_locus: call 0"]
    M -->|Resolved candidate sequence| O["Update call and evidence state"]
    M -->|Unresolved or case limit reached| P["Retain evidence and mark unresolved or deferred"]
    K -->|Initial negative with no evidence| Q["not_detected: call 0; absence not proven"]
```

The flowchart summarizes routing. The evidence table below defines the primary
states; arbitration can refine a provisional decision.

| Evidence state | Default interpretation | Binary treatment |
| --- | --- | --- |
| `confirmed_present` | Breadth ≥95%, identity ≥95%, and sufficient depth. | `1` |
| `divergent_variant` | Breadth ≥90%, identity ≥90% but <95%, and sufficient depth. | `1`, flagged as divergent |
| `possible_truncation` | Breadth ≥70% but <95%, identity ≥95%, and sufficient depth. | Initial positive remains provisionally `1`; otherwise unresolved |
| `partial_homolog` | Partial coverage or weak sequence similarity does not support an intact gene. | `0`, with evidence retained |
| `ambiguous_multimap` | Reads support a family but cannot resolve the exact cluster. | Unresolved; no automatic presence for every homolog |
| `not_detected` | No meaningful read evidence. | Initial negative stays `0`; initial positive requires arbitration |
| `insufficient_evidence` | Depth or reconstructed identity is insufficient for a decision. | Unresolved; preserve the initial binary call unless resolved |
| `confirmed_absent_locus` | A local reconstruction supports a deletion junction spanning both flanks. | `0`, with physical absence evidence |

**Zero unique mappings do not prove biological absence.** Unresolved validated
calls can be blank in the evidence table. When no resolved replacement exists,
the binary matrix preserves the initial call; consult the evidence table when
interpreting that value. `partial_homolog` and `confirmed_absent_locus` both permit
an intact-gene call of zero but describe different biology.

Identity is measured from reconstructed sequence rather than average read
alignment identity. Own-locus identity uses the sample CDS as its reference;
recovery uses the pangenome reference. Normalized depth is gene depth divided by
a representative sample depth estimate and is supporting evidence, not a universal
presence threshold.

## Installation

Requirements: Linux, Git, and a working Miniforge/Mambaforge installation providing
`mamba`. Slurm is required only for the Slurm execution profile.

```bash
git clone https://github.com/AndriyPlakhotnyk/CleanGene.git
cd CleanGene
bash scripts/install_or_update.sh
conda activate cleangene
```

The installer creates or updates the bioinformatics environments, installs
CleanGene from the checkout, and verifies its tools. CheckM2 runs through an
automatically managed companion environment. Installing the Python package alone
with `pip` does not install the external bioinformatics tools.

To update your current branch:

```bash
git pull --ff-only
bash scripts/install_or_update.sh
conda activate cleangene
```

## Quick start

### 1. Configure the run

The installer creates `config/cleangene.arc.local.env` from the supplied Slurm
template without overwriting an existing local file. Set `SLURM_ACCOUNT` and
`SLURM_PARTITION` for your cluster, and adjust resource requests as needed.

```bash
cleangene doctor --config config/cleangene.arc.local.env
```

Managed databases are downloaded when needed and reused. Set
`CLEANGENE_DATABASE_ROOT` to place them on a shared filesystem; use `KRAKEN2_DB`
or `CHECKM2_DB` to select an existing database explicitly.

### 2. Prepare a manifest

Use a tab-separated file with one row per isolate. The example below contains
paired FASTQ paths and an explicit analysis group:

```tsv
isolate_id	group_id	R1	R2
isolate_01	species_A	/data/isolate_01_R1.fastq.gz	/data/isolate_01_R2.fastq.gz
isolate_02	species_A	/data/isolate_02_R1.fastq.gz	/data/isolate_02_R2.fastq.gz
```

Use groups of biologically comparable isolates. Pangenome analysis requires at
least two retained isolates per group.

| Column | Purpose |
| --- | --- |
| `isolate_id` | Unique sample identifier. |
| `R1`, `R2` | Paired FASTQ inputs; compressed FASTQs are supported. |
| `raw_bam` | Alternative to `R1`/`R2`: a paired, unmapped sequencing BAM containing both mates. |
| `group_id` | Explicit pangenome group. |
| `organism` | Defines the group when `group_id` is absent. Without either field, Kraken2 determines grouping. |
| `assembly`, `gff` | Reusable assembly and annotation artifacts. |
| `pangenome_dir` | Existing Panaroo output to use instead of generating a new pangenome. |

Provide either paired FASTQs or a uBAM for each isolate. Optional reuse fields
include `reads_processed`, `fastp_json`, `kraken_report`, `protein_fasta`, and
`checkm2_report`.

### 3. Launch

```bash
cleangene run \
  --profile slurm \
  --manifest input/cohort.manifest.tsv \
  --analysis-root /data/cleangene-analysis \
  --config config/cleangene.arc.local.env \
  --assembler spades
```

This submits a controller job and returns to the shell. Add `--ignore-checkm2`
to omit CheckM2 completeness and contamination assessment. Other QC criteria
remain active.

For a small local analysis, use `--profile local`. For a Slurm submission preview,
add `--dry-run`; it creates run metadata and prints the controller submission
command without submitting jobs. **Dry-run is a Slurm option; do not use it to
preview local execution.**

## Command-line arguments

Run `cleangene run --help` for the complete CLI reference. Configuration values
apply unless overridden by a command-line option.

| Argument | Values / default | Purpose |
| --- | --- | --- |
| `--manifest` | TSV path | Sample inputs; required for a new run. |
| `--analysis-root` | Directory | Parent directory for `runs/<run-id>/`. Required by `run`. |
| `--config` | Environment-style file | QC, database, and execution settings. |
| `--profile` | `slurm` (default), `local` | Execution backend. |
| `--assembler` | `shovill` (built-in default), `spades`, `off` | Assembly strategy; `off` skips assembly and annotation. |
| `--ignore-checkm2` | Flag | Skip CheckM2 assessment during the run. |
| `--skip-trim` | Flag | Bypass fastp trimming. |
| `--compress-assembly-outputs` | `off`, `intermediates`, `all` | Assembly storage policy; built-in default is `intermediates`. |
| `--compress-annotation-outputs` | `off`, `nonessential` | Annotation storage policy; built-in default is `nonessential`. |
| `--cleanup-trimmed-fastq` | Flag | Enable final cleanup of retained trimmed FASTQs. |
| `--run-id` | Identifier | Set a run name instead of the generated timestamp. |
| `--resume` | Existing run ID | Resume a run under the analysis root. |
| `--dry-run` | Flag, Slurm profile | Print submission command without submitting. |
| `--cancel-active` | Flag, resume | Cancel active jobs associated with the run before resubmission. |

Assembly intermediates and nonessential annotation outputs are compressed by
default; the primary assembly and GFF remain available. Explicit `off` settings
in existing configs still override these defaults. To adopt the new policy in
an existing local config, set `COMPRESS_ASSEMBLY_OUTPUTS="intermediates"` and
`COMPRESS_ANNOTATION_OUTPUTS="nonessential"`.

Direct SPAdes mode uses original paired reads with `--only-assembler`. Choose
Shovill when you want its assembly preparation workflow.

### Validation settings

Set these in your config file. Breadth and identity thresholds use fractions
between zero and one, not percentages.

| Setting | Built-in default | Purpose |
| --- | --- | --- |
| `VALIDATION_SCOPE` | `all` | Validate all genes; `accessory` and `differential` restrict selection. |
| `READ_VALIDATION_MIN_BREADTH` | `0.95` | Confirmed-presence breadth. |
| `READ_VALIDATION_MIN_IDENTITY` | `0.95` | Confirmed-presence sequence identity. |
| `READ_VALIDATION_MIN_MEAN_DEPTH` | `5` | Minimum mean depth for sequence-based calls. |
| `READ_VALIDATION_MIN_MAPQ` | `20` | High-confidence mapping threshold. |
| `READ_VALIDATION_TRUNCATION_MIN_BREADTH` | `0.70` | Lower breadth threshold for possible truncation. |
| `READ_VALIDATION_DIVERGENT_MIN_BREADTH` | `0.90` | Minimum divergent-variant breadth. |
| `READ_VALIDATION_DIVERGENT_MIN_IDENTITY` | `0.90` | Minimum divergent-variant identity. |
| `READ_VALIDATION_ARBITRATION_MAX_CASES` | `20` | Maximum reconstruction cases per isolate. |
| `READ_VALIDATION_ARBITRATION_MAX_READS` | `100000` | Maximum recruited read names per reconstruction. |
| `READ_VALIDATION_ARBITRATION_MEMORY_GB` | `12` | SPAdes memory cap for reconstruction. |
| `READ_VALIDATION_FLANK_LENGTH` | `500` | Flanking bases used for locus investigation. |
| `READ_VALIDATION_DELETION_MIN_IDENTITY` | `0.95` | Minimum junction sequence identity. |
| `READ_VALIDATION_DELETION_MIN_ANCHOR` | `50` | Contiguous aligned bases on each side of the junction. |

Use `SLURM_PREPROCESS_MAX_INFLIGHT`, `SLURM_VALIDATION_MAX_INFLIGHT`, and
`SLURM_ARBITRATION_MAX_INFLIGHT` to control concurrency. Stage-specific CPU,
memory, and time requests are available in the
[Slurm configuration template](config/cleangene.arc.env). The
[example config](config/cleangene.example.env) and
[default settings](src/cleangene/defaults.py) list additional controls, including
QC thresholds and database management.

## Outputs

Each run is stored beneath `<analysis-root>/runs/<run-id>/`.

| Location within a run | Contents |
| --- | --- |
| `results/groups/<group>/cleaned_pangenome.tsv` | Final binary gene presence/absence matrix. |
| `results/groups/<group>/03_read_validation/gene_call_evidence.long.tsv` | Initial/final calls, evidence states, resolution, and arbitration details. |
| `results/groups/<group>/03_read_validation/read_validation_metrics.tsv` | Coverage, depth, normalized depth, identity, reconstruction lengths, mappings, and CDS coordinates. |
| `results/sample_data/<isolate>/` | Per-isolate QC and processing artifacts. |
| `results/cohort/` | Cohort QC and summaries. |
| `provenance/` | Manifest, resolved configuration, and runtime metadata. |
| `logs/slurm/` | Controller and stage job logs. |

### Validation summaries

Each group's `03_read_validation/` directory also contains:

| File | Contents |
| --- | --- |
| `gene_call_summary.tsv` | Group totals: all calls, kept, changed, 0→1, 1→0, initial/final presence counts, and percentages. |
| `gene_call_summary.per_isolate.tsv` | The same counts and percentages for every retained isolate. |
| `summary_statistics.txt` | Reconciled final-matrix statistics, including core, soft-core, shell, cloud, and absent-in-all clusters. |
| `decision_reason_upset.png`, `.svg` | UpSet-style plot: decision reasons across columns, change direction as color, and deciding metrics in the bottom matrix. |
| `decision_reason_upset.tsv` | Exact counts, metric membership, and provenance used to draw the plot. |

“Kept” means an unchanged call, including both 0→0 and 1→1. An isolate's percentage
denominator is **all gene clusters in its group's matrix**, including untested
calls carried forward. Cohort totals count gene-isolate calls across groups;
percentages are weighted by those counts. Cohort tables and the combined plot are
written to `results/cohort/`, including `gene_call_summary.per_group.tsv`.

The text report is published only after initial/final matrices and evidence
agree, kept + changed equals total calls, and both change directions reconcile.
Its prevalence bins are recalculated from final calls; provisional calls remain
part of the binary matrix. Cloud counts include clusters absent from all retained
isolates, which are also listed separately.

Plot dots mean **metrics used by the deciding rule**, not thresholds passed.
New evidence records `decision_metrics`; older evidence uses labeled state-based
inference. Reasons with different metric sets appear in separate columns. Reports
with more than 30 columns have additional numbered pages; no decisions are omitted.
A run with no changed calls receives an explicit “No gene calls changed” plot.
Resume backfills missing reports using existing validation evidence, without
rerunning mapping or arbitration solely for reporting.

Keep the evidence tables with the binary matrix. They identify provisional calls,
family-only evidence, partial homologs, and supported deletions that a binary value
cannot express.

## Resume and monitor

```bash
cleangene resume \
  --run-dir /data/cleangene-analysis/runs/<run-id> \
  --config config/cleangene.arc.local.env
```

Completed work is reused where its completion checks pass. Resume reconciles
preprocessing outputs and reruns incomplete or invalidated stages. For local
resumption, use `cleangene run --profile local --analysis-root ... --resume <run-id>`.
Use a new run when changing biological thresholds.

Monitor Slurm jobs with `squeue -j <job-id>` and inspect the run's controller and
stage logs. To audit missing preprocessing completion markers:

```bash
cleangene reconcile-preprocess --run-dir /data/cleangene-analysis/runs/<run-id>
```

## Downstream analysis

`cleangene utils` provides sample selection, differential gene analysis, operon
analysis, read-backed variants, iTOL exports, and post-hoc CheckM2 assessment.
See the [utilities guide](docs/UTILS.md) for commands and output descriptions.

## Interpretation and scope

Local reconstruction is limited to discordant cases. Exact allele resolution
among highly similar homologs can remain unresolved, and automatic transfer of
flanking loci from other isolates is not currently implemented. ORF checks are
basic; normalized depth uses an assembly-based chromosomal proxy. Cohort-scale
performance depends on reference complexity, coverage, and available resources.

## Testing

With the CleanGene environment active:

```bash
bash tests/run_tests.sh
```

The suite covers classification, local and Slurm orchestration, resume, QC, and
storage behavior. Synthetic-read integration tests exercise mapping, recovery,
local assembly, and deletion-junction detection when their tools are available.

## License

See [LICENSE](LICENSE).
