# CleanGene

CleanGene is a Slurm-native bacterial isolate workflow for read QC, assembly,
annotation, pangenome construction, read-backed gene validation, and cohort QC.

## Installation and update

The supported ARC installation is:

```bash
git clone https://github.com/AndriyPlakhotnyk/CleanGene.git
cd CleanGene
bash scripts/install_or_update.sh
conda activate cleangene
```

For an existing checkout:

```bash
git pull
bash scripts/install_or_update.sh
conda activate cleangene
```

The script creates missing environments, updates existing environments in
place, installs the current checkout, preserves an existing local config, and
runs deployment checks. To remove and rebuild both managed environments:

```bash
bash scripts/install_or_update.sh --recreate
conda activate cleangene
```

CheckM2 is isolated in a companion environment because its TensorFlow and
DIAMOND dependencies conflict with parts of the primary bioinformatics
toolchain. CleanGene creates and calls that environment automatically. Do not
activate it or configure its executable during normal use.

`pip install .` and `pip install -e .` install only the CleanGene Python
package. They do not install Shovill, SPAdes, Prokka, Panaroo, Kraken2, CheckM2,
or the other required bioinformatics tools. Use the setup script for a complete
installation.

## ARC configuration

The setup script creates `config/cleangene.arc.local.env` from the tracked
`config/cleangene.arc.env` template only when the local file is absent. It never
overwrites an existing local config.

Edit site-specific values in the ignored local file, then verify them:

```bash
vi config/cleangene.arc.local.env
cleangene doctor --config config/cleangene.arc.local.env
```

Leave `CHECKM2_EXECUTABLE`, `CHECKM2_DB`, `CHECKM2_DATABASE_ROOT`, and
`KRAKEN2_DB` blank for automatic software and database management. Set
`CLEANGENE_DATABASE_ROOT` only when managed databases should live on a shared
filesystem outside the checkout. Keep account names, private paths, and other
site-specific values out of the tracked template.

## Manifest

A minimal paired-FASTQ manifest is tab separated:

```text
isolate_id  group_id  R1                    R2
ERR001      sp00001   /data/r1.fastq.gz     /data/r2.fastq.gz
```

Raw sequencing BAM containers are accepted instead of FASTQ pairs:

```text
isolate_id  group_id  raw_bam
ERR001      sp00001   /data/raw_reads.bam
```

Grouping behavior:

- A supplied `group_id` is used directly.
- With no `group_id`, a supplied `organism` defines the group.
- With neither value, Kraken2 identifies the top species and defines the group.
- A supplied `pangenome_dir` skips pangenome generation and validates against
  that existing Panaroo output.

Reusable artifact columns include `reads_processed`,
`read_processing_pipeline`, `read_processing_version`, `read_qc_tsv`,
`fastp_json`, `kraken_report`, `assembly`, `gff`, `protein_fasta`,
`checkm2_report`, `expected_genome_size`, `prodigal_training_file`, and
`pangenome_dir`.

## Dry run

```bash
cleangene run \
    --manifest input/manifest.tsv \
    --analysis-root /path/to/cleangene-output \
    --config config/cleangene.arc.local.env \
    --profile slurm \
    --dry-run
```

The dry run prints the `sbatch` commands without submitting jobs.

## Run

```bash
cleangene run \
    --manifest input/manifest.tsv \
    --analysis-root /path/to/cleangene-output \
    --config config/cleangene.arc.local.env
```

CleanGene resolves and records runtime executables and database paths before
or during their dedicated setup stages. Missing managed databases are prepared
once and reused by every isolate and later run.

## Resume and monitoring

```bash
cleangene resume \
    --run-dir /path/to/cleangene-output/runs/<run-id> \
    --config config/cleangene.arc.local.env
```

Alternatively, use `cleangene run --analysis-root ... --resume <run-id>`.
Completed state markers, resolved database paths, and the resolved CheckM2
executable are retained.

Monitor jobs with `squeue -j <job-id>`. Controller and stage logs are under
`<run>/logs/slurm/`; per-isolate Slurm logs are grouped in
`logs/slurm/preprocess/` and `logs/slurm/validate/`. CheckM2 stdout, stderr, and
elapsed prediction time are recorded with each isolate's preprocess logs.

On resume, CleanGene reconciles missing preprocess markers against terminal
`results/sample_data/*/qc.tsv` outputs before submitting work. Ambiguous,
malformed, or incomplete evidence fails closed. Audit or repair markers with:

```bash
cleangene reconcile-preprocess --run-dir /path/to/run
cleangene reconcile-preprocess --run-dir /path/to/run --apply --require-all
```

Before downstream pangenome work starts, exclude isolates without changing
manifest row indices:

```bash
cleangene exclude --run-dir /path/to/run --samples-file unwanted_samples.txt
```

## Scientific and QC details

CleanGene classifies each isolate as `PASS`, `WARNING`, or `FAIL`. Boundary
values belong to the less severe state.

| Metric | PASS | WARNING | FAIL/exclude |
|---|---|---|---|
| Expected organism | Kraken top species matches | Classification unavailable | Top species differs |
| Kraken foreign-species contamination | `<=5%` | None | `>5%` |
| CheckM2 completeness | `>=90%` | `>=80%` and `<90%` | `<80%` |
| CheckM2 contamination | `<=5%` | `>5%` and `<=10%` | `>10%` |
| Assembly contigs | `<=300` | `301-1000` | `>1000` |
| Assembly N50 | `>=25000 bp` | `5000-24999 bp` | `<5000 bp` |
| Sequencing coverage | `>=20x` | `>=10x` and `<20x` | `<10x` |
| Post-processing read length | `>=120 bp` | `<120 bp` | No default hard fail |
| Post-processing mean base quality | `>=Q30` | `<Q30` | No default hard fail |
| Internal assembly and Prokka GFF | Both produced | Not evaluated for an external pangenome | Missing or failed |
| Explicit user exclusion | - | - | `exclude=true` or `user_excluded=true` |

The canonical thresholds are:

```text
QC_MAX_CONTIGS_PASS=300
QC_MAX_CONTIGS_FAIL=1000
QC_MIN_N50_PASS=25000
QC_MIN_N50_FAIL=5000
QC_MIN_COVERAGE_PASS=20
QC_MIN_COVERAGE_FAIL=10
QC_MIN_READ_LENGTH_PASS=120
QC_MIN_READ_LENGTH_FAIL=
QC_MIN_MEAN_BASE_QUALITY_PASS=30
QC_MIN_MEAN_BASE_QUALITY_FAIL=
QC_MIN_COMPLETENESS_PASS=90
QC_MIN_COMPLETENESS_FAIL=80
QC_MAX_CHECKM2_CONTAMINATION_PASS=5
QC_MAX_CHECKM2_CONTAMINATION_FAIL=10
QC_MAX_KRAKEN_CONTAMINATION_FAIL=5
QC_PROFILE_FILE=
```

Blank read-length and mean-quality fail thresholds make those criteria
warning-only. Numeric fail thresholds enable normal PASS, WARNING, and FAIL
bands. Every global threshold has a matching lowercase manifest override.
Blank manifest cells inherit lower-priority values.

`QC_PROFILE_FILE` accepts a TSV with `scope_type`, `scope_value`, and lowercase
threshold columns. Precedence is per-isolate manifest override, matching group
profile, matching organism profile, global config, then built-in default.
Resolved thresholds are copied into run provenance.

`trimmed_read_length` is the smaller mean R1/R2 length. `mean_base_quality` is
the weighted mean Phred score across all bases. `sequencing_coverage` is total
post-processing read bases divided by assembly length. CheckM2 values come from
`quality_report.tsv`. WARNING isolates continue downstream; FAIL isolates are
excluded.

## Advanced configuration

### Runtime architecture

```text
cleangene environment
    ├── Shovill / SPAdes
    ├── Prokka
    ├── Panaroo
    ├── BWA / samtools / bcftools
    ├── Kraken2
    └── CleanGene
             |
             └── automatically calls
                 cleangene-checkm2/checkm2
                         |
                         └── shared managed CheckM2 DB
```

The CheckM2 executable is resolved before controller submission. Resolution
honors an explicit `CHECKM2_EXECUTABLE`, then `PATH`, an executable beside the
active Python, and finally the sibling `cleangene-checkm2` environment. The
absolute executable and version are stored in `provenance/resolved_config.json`
and used for `--version`, database download, and prediction.

With normal settings, `checkm2_db_setup` reuses
`CheckM2_database/uniref100.KO.1.dmnd` below the managed root or downloads it
once under an inter-process lock. Download uses the resolved companion
executable with `--no_write_json_db`; prediction always supplies
`--database_path` and `--remove-intermediates`. An invalid explicit
`CHECKM2_DB` fails closed.

`CHECKM2_CPUS`, `CHECKM2_MEM`, and `CHECKM2_TIME` size database setup.
Per-isolate prediction runs inside preprocess and uses `CPUS`. A supplied
`checkm2_report` remains reusable. `CHECKM2_MODE=off` is supported but produces
a QC warning because completeness and contamination were not evaluated.

### Managed Kraken2 databases

With `KRAKEN2_DB=""`, CleanGene manages `kraken2_<size>` below
`CLEANGENE_DATABASE_ROOT` or the checkout's `databases/` directory. Supported
sizes are `standard-8`, `standard-16`, and `standard`. A valid database contains
non-empty `hash.k2d`, `opts.k2d`, and `taxo.k2d`. Setup is locked and shared
across runs.

`KRAKEN2_DATABASE_ROOT` overrides only the Kraken2 parent. `KRAKEN2_DB` is an
exact custom database override. `CHECKM2_DATABASE_ROOT` and `CHECKM2_DB` provide
the corresponding CheckM2-specific overrides.

### Preprocessing modes

`--assembler shovill` is the default. `--assembler spades` runs direct SPAdes
with original untrimmed paired reads and `--only-assembler`. `--assembler off`
or legacy `--skip_shovill` performs read/Kraken QC without assembly or
annotation. `--skip_trim` bypasses fastp while retaining assembly.

Storage controls:

```text
COMPRESS_ASSEMBLY_OUTPUTS=off|intermediates|all
COMPRESS_ANNOTATION_OUTPUTS=off|nonessential
CLEANUP_TRIMMED_FASTQ=false|true
```

Cleanup can also be run after completion:

```bash
cleangene cleanup --run-dir /path/to/run --dry-run
cleangene cleanup --run-dir /path/to/run
```

Original input files are never modified.

### Scale and outputs

The ARC template maintains rolling, capacity-aware arrays and prioritizes known
groups from smallest to largest. New runs use an indexed isolate task store so
workers seek directly to their manifest records. Job-count and optional CPU
headroom controls prevent oversubmission.

Final group matrices are written to
`results/groups/<group>/cleaned_pangenome.tsv`. Sample-specific outputs live in
`results/sample_data/<isolate>/`. `results/organisms/<organism>/<isolate>` is a
storage-free symlink index into sample data. Cohort summaries are under
`results/cohort/`.

### Local debugging and utilities

Small local tests can use `--profile local`. This does not replace the supported
ARC installation and is not intended for large cohorts.

Completed runs provide Slurm-native sample queries, differential gene tests,
operon typing, read-backed variant analysis, and iTOL datasets through
`cleangene utils`. See [docs/UTILS.md](docs/UTILS.md).
