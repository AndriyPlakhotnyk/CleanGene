
## Quick Start on a SLURM cluster

1. **Install Conda/Rocks**   
   ```bash
   mamba env create -f environment.yml
   mamba env create -f environment.checkm2.yml
   conda activate cleangene
   pip install -e .
   ```

   CheckM2 uses a companion environment because its TensorFlow and DIAMOND
   dependencies conflict with parts of the Prokka/Shovill toolchain. CleanGene
   automatically discovers `cleangene-checkm2` beside the active `cleangene`
   environment; users do not need to activate it or configure an executable
   path. After pulling environment changes, update both environments with:

   ```bash
   mamba env update --name cleangene --file environment.yml --prune
   mamba env update --name cleangene-checkm2 --file environment.checkm2.yml --prune
   conda activate cleangene
   pip install -e .
   cleangene doctor --config config/cleangene.arc.local.env
   ```

2. **Prepare configuration**  
   Copy the example and fill your paths. Keep machine-specific values in an
   ignored local config when possible, especially accounts, private database
   roots, and exact executable paths.

   ```bash
   cp config/cleangene.example.env config/cleangene.env
   vi config/cleangene.env
   ```

   For ARC or another shared cluster config, prefer:

   ```bash
   cp config/cleangene.arc.env config/cleangene.arc.local.env
   vi config/cleangene.arc.local.env
   cleangene doctor --config config/cleangene.arc.local.env
   ```

3. **Create a manifest**  
   The `input/manifest.tsv` should list isolates and FASTQ pairs, e.g.:  

   ```text
   isolate_id  group_id  R1         R2
   ERR001      sp00001   /data/r1.fastq.gz  /data/r2.fastq.gz
   ```

   Raw sequencing BAM containers are also accepted instead of FASTQ pairs:

   ```text
   isolate_id  group_id  raw_bam
   ERR001      sp00001   /data/raw_reads.bam
   ```

   Grouping rules:

   - `group_id` supplied: use that pangenome group.
   - no `group_id`, `organism` supplied: group by organism.
   - neither supplied: run Kraken2 and group by inferred top species.
   - `pangenome_dir` supplied: skip pangenome generation and run BWA validation against that Panaroo output.

4. **Run a test dry‑run**

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.env \
       --profile slurm --dry-run
   ```

   The command prints the `sbatch` commands that would be submitted. Verify they appear reasonable and your cluster accepts them.

5. **Submit the real workflow**

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.env
   ```

6. **Monitoring**  
   Job status with `squeue -j <JOB_ID>`. Controller and downstream logs are in
   `<run>/logs/slurm/`; individual preprocess and read-validation logs are
   grouped under `<run>/logs/slurm/preprocess/` and
   `<run>/logs/slurm/validate/`.

7. **Large scale (≥30 k isolates)**  
   Use `config/cleangene.arc.env` on ARC. It leaves `SLURM_PARTITION` empty and maintains a rolling, capacity-aware preprocessing window of at most 400 submitted tasks. Up to eight array batches may be active, but every submission still respects total user occupancy and QOS headroom. Known organism groups are prioritized smallest first and can enter Panaroo/validation while later groups are still preprocessing. ARC config sets `TAXONOMY_MODE=auto` and `KRAKEN2_BUNDLE_SIZE=8`.

   New runs write an indexed isolate task store at
   `state/tasks/isolate_tasks.jsonl` plus `state/tasks/isolate_tasks.idx`.
   Each isolate worker seeks directly to its own JSON record instead of
   scanning the full manifest or QC-threshold table. Older runs are migrated by
   the SLURM controller during resume without deleting prior evidence.

   Optional manifest artifact columns are preserved for reuse and provenance:
   `reads_processed`, `read_processing_pipeline`, `read_processing_version`,
   `read_qc_tsv`, `fastp_json`, `kraken_report`, `assembly`, `gff`,
   `protein_fasta`, `checkm2_report`, `expected_genome_size`, and
   `prodigal_training_file`. Path validation happens in worker jobs, not during
   login-node submission.

   CleanGene automatically manages Kraken2 databases when Kraken2 taxonomy is
   needed. With `KRAKEN2_DB=""` and
   `KRAKEN2_DATABASE_SIZE="standard-8"`, CleanGene checks:

   ```text
   <CleanGene>/databases/kraken2_standard-8/
   ```

   A valid database must contain non-empty `hash.k2d`, `opts.k2d`, and
   `taxo.k2d`. If the selected shared database exists, it is reused. If it is
   absent and `KRAKEN2_AUTO_DOWNLOAD=true`, CleanGene downloads it once into
   that shared directory; later runs reuse the same path. Supported sizes are
   `standard-8`, `standard-16`, and `standard`; aliases such as `8`, `16`,
   `standard_8`, `standard_16`, `full`, and `full-standard` normalize to those
   canonical names.

   Set `CLEANGENE_DATABASE_ROOT=/path/to/shared/data/CleanGene` on HPC systems
   that keep large databases on a dedicated shared filesystem. CleanGene will
   manage selected Kraken2 databases below that parent, for example
   `/path/to/shared/data/CleanGene/kraken2_standard-8`. `KRAKEN2_DATABASE_ROOT`
   can override only the Kraken2 parent, and `KRAKEN2_DB` remains an exact
   database-directory override for intentionally using a custom Kraken2
   database. If CleanGene is installed in a layout where the source/project root
   cannot be identified, the automatic root falls back to
   `~/.cache/cleangene/databases`; set `CLEANGENE_DATABASE_ROOT` to make the
   shared location explicit.

8. **FASTQ QC only / skip trimming**

   To run only read-level QC for supplied FASTQs and avoid assembly storage,
   skip Shovill. CleanGene records isolate QC and stores stable run-local
   symlinks to the original FASTQ files for resume and utility lookups.

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.arc.env \
       --skip_shovill \
       --skip_trim
   ```

   Equivalent config options are `SKIP_SHOVILL=true` and `SKIP_TRIM=true`.
   `--skip_trim` is also useful on full runs when you want to skip fastp
   adapter trimming but still assemble with Shovill from the original reads.
   With `READ_TRIMMING_MODE=auto`, manifest `reads_processed=true` skips fastp
   explicitly; `READ_TRIMMING_MODE=always` still runs fastp. Isolate `qc.tsv`
   includes `read_processing_decision` explaining why trimming ran or was
   skipped.

   `--skip_shovill` is the legacy QC-only option and skips all assembly,
   including SPAdes. To bypass Shovill but still assemble the original,
   untrimmed paired FASTQs directly with SPAdes, use:

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.arc.env \
       --assembler spades
   ```

   `--assembler shovill` is the default full Shovill workflow;
   `--assembler off` is equivalent to `--skip_shovill`. Direct SPAdes implies
   trimming off, invokes `spades.py --only-assembler` to skip SPAdes read-error
   correction, and stores run-local symlinks to the original FASTQs.

   To reduce Shovill assembly storage without skipping assembly, use:

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.arc.env \
       --compress_assembly_outputs intermediates
   ```

   `COMPRESS_ASSEMBLY_OUTPUTS=intermediates` gzips Shovill/SPAdes leftover
   FASTA/GFA files such as `spades.fasta` and `spades.gfa`, while keeping the
   final `contigs.fa`/`contigs.fasta` plain. `COMPRESS_ASSEMBLY_OUTPUTS=all`
   also gzips the final contig file and records the `.gz` path in isolate QC.

   To compress nonessential Prokka annotation outputs while keeping `.gff`
   available for Panaroo, resume checks, and utilities, set:

   ```bash
   COMPRESS_ANNOTATION_OUTPUTS=nonessential
   ```

   or pass `--compress_annotation_outputs nonessential` to `run`/`resume`.

9. **Resume an unfinished run**

   ```bash
   cleangene run \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.env \
       --resume <run-id>
   ```

   `cleangene resume --run-dir /path/to/run --config config/cleangene.arc.env`
   refreshes execution settings for an existing run while retaining completed
   task markers, resolved database paths, and the resolved CheckM2 executable
   when present.

   For a storage-saving resume after all orphaned SLURM tasks have stopped, use:

   ```bash
   cleangene resume \
       --run-dir /path/to/run \
       --config config/cleangene.arc.env \
       --cleanup_trimmed_fastq \
       --compress_assembly_outputs intermediates \
       --compress_annotation_outputs nonessential
   ```

   The final summary job also sweeps outputs from preprocess tasks that finished
   before the resumed controller was submitted. It compresses safe run-local
   assembly intermediates and nonessential annotation files, then replaces
   CleanGene's trimmed FASTQs with symlinks to the original manifest FASTQs.
   The original inputs are opened only as symlink targets and are never modified.

   On resume, the SLURM controller reconciles missing preprocess markers against
   existing terminal isolate QC outputs before submitting new preprocess arrays.
   It scans `results/sample_data/*/qc.tsv` and legacy
   `results/groups/*/01_isolates/*/qc.tsv` outputs once, writes
   `state/preprocess_reconciliation.tsv`, and recreates missing
   `state/preprocess/<isolate>.done.json` markers only when the QC row is a
   valid terminal result. A retained internal PASS/WARNING isolate must have
   readable nonzero assembly and GFF files; terminal FAIL/excluded isolates,
   external-pangenome isolates, and `ASSEMBLER=off` isolates do not require
   assembly/GFF for marker recovery. The current group path is checked first;
   if it has no QC candidate, historical groups such as `__kraken_pending__`
   are accepted only when exactly one candidate validates. Ambiguous, malformed,
   multirow, wrong-isolate, stale-marker, or legacy QC output fails closed before
   any new `sbatch` submission. Suspicious existing markers are quarantined as
   `*.stale.<timestamp>` instead of being silently overwritten.

   Example controller summary:

   ```text
   step=preprocess_reconciliation | total=30000 | marker_complete=26000 | output_recovered=3465 | active=20 | incomplete=515 | inconsistent=0
   ```

   To deliberately reprocess inconsistent isolates without deleting their old
   outputs, set `RESUME_REPROCESS_INCONSISTENT_PREPROCESS=true`. The default is
   `false` so questionable evidence fails closed.

   You can audit or repair marker state manually without submitting the full
   controller:

   ```bash
   cleangene reconcile-preprocess --run-dir /path/to/run
   cleangene reconcile-preprocess --run-dir /path/to/run --apply --require-all
   ```

   The first command is a dry run. The second recreates only markers backed by
   valid outputs and exits nonzero if any isolate remains incomplete or
   inconsistent. `--compress-safe` may be added with `--apply` to run the same
   safe assembly/annotation compression sweep.

   ARC users who also need a CPU-based submit cap can set
   `SLURM_USER_CPU_LIMIT` and `SLURM_CPU_HEADROOM`; both default to `0`, which
   preserves the existing job-count-only behavior.

   Do not delete rows from a run's `provenance/manifest.tsv`: array task indices
   and completed state markers must remain stable. Before any downstream
   pangenome stage starts, exclude samples safely with:

   ```bash
   cleangene exclude --run-dir /path/to/run --samples-file unwanted_samples.txt
   ```

   The file contains one isolate ID per line. This marks rows as
   `user_excluded`, retains provenance and task indices, skips unfinished work
   for those isolates, and removes already-preprocessed isolates from Panaroo
   and read validation. If downstream work has already started, CleanGene
   refuses the change and a new filtered run is required.

   Each completed organism/group publishes its final read-validated presence/absence matrix at `results/groups/<group>/cleaned_pangenome.tsv`. Panaroo intermediates and the detailed BWA validation evidence remain in the numbered subdirectories.

   After organism resolution, CleanGene also creates a storage-free isolate
   index at `results/organisms/<organism>/<isolate_id>`. Each isolate entry is a
   relative symlink to its complete directory under `results/sample_data`, exposing
   its annotation, assembly, reads, QC, and per-isolate logs without copying
   them. `results/cohort/organism_isolate_index.tsv` records the display names,
   symlink paths, and canonical targets. The index is refreshed during final
   summary and is safe to regenerate when resuming a run.

10. **Local debugging** (small tests)

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ./local_run \
       --profile local
   ```

11. **Downstream analyses**

   Completed runs can submit SLURM-native sample queries, differential-gene tests, operon typing, read-backed variant analysis, and iTOL annotation datasets:

   ```bash
   cleangene utils get-samples \
       --run-dir ~/cleangene-output/runs/<run-id> \
       --organism "Species name" --genes geneA geneB
   ```

   See [docs/UTILS.md](docs/UTILS.md) for command forms, manifest columns, outputs, and resource controls.

12. **Reclaim trimmed-read storage after a completed run**

   Preview and then replace generated trimmed FASTQ pairs with symlinks to the
   original manifest FASTQs. The stable paths recorded in `qc.tsv` continue to
   work for read-backed `cleangene utils` commands.

   ```bash
   cleangene cleanup --run-dir /path/to/run --dry-run
   cleangene cleanup --run-dir /path/to/run
   ```

   For automatic cleanup at the end of future runs, pass
   `--cleanup_trimmed_fastq` or set `CLEANUP_TRIMMED_FASTQ=true` in the run
   configuration. Cleanup occurs in the
   final summary job, after all core read consumers finish. FASTQs extracted
   from raw BAM inputs are retained because they have no original FASTQ target.
   Details and reclaimed byte counts are written to
   `results/cohort/fastq_cleanup.tsv`.

13. **Isolate QC classification**

   CleanGene classifies every isolate as `PASS`, `WARNING`, or `FAIL` using the
   following exact boundaries. Boundary values belong to the less severe
   state.

   | Metric | PASS | WARNING | FAIL/exclude |
   |---|---|---|---|
   | Expected organism | Kraken top species matches expected organism | Expected organism supplied but classification unavailable | Normalized top species differs from expected organism |
   | Kraken foreign-species contamination | `<=5%` | None | `>5%` |
   | CheckM2 completeness | `>=90%` | `>=80%` and `<90%` | `<80%` |
   | CheckM2 contamination | `<=5%` | `>5%` and `<=10%` | `>10%` |
   | Assembly contigs | `<=300` | `301–1000` | `>1000` |
   | Assembly N50 | `>=25000 bp` | `5000–24999 bp` | `<5000 bp` |
   | Sequencing coverage | `>=20x` | `>=10x` and `<20x` | `<10x` |
   | Post-processing read length | `>=120 bp` | `<120 bp` | No default hard fail |
   | Post-processing mean base quality | `>=Q30` | `<Q30` | No default hard fail |
   | Assembly and Prokka GFF | Both produced when CleanGene builds the pangenome | Not evaluated for an explicitly supplied external pangenome | Missing or failed during an internal pangenome run |
   | Explicit user exclusion | — | — | `exclude=true` or `user_excluded=true` |

   Organism comparison is case-insensitive and collapses whitespace. If no
   expected `organism` was supplied, taxonomic agreement is not evaluated.
   `trimmed_read_length` is the smaller of the mean R1 and R2 lengths for the
   actual post-processing files passed to assembly. `mean_base_quality` is the
   weighted mean Phred score across every base in both files; it is not the
   percentage of Q30 bases. `sequencing_coverage` is total post-processing R1
   and R2 bases divided by assembly length. Values are not rounded before
   classification. Valid `fastp.json` length/base statistics are reused, and
   the FASTQs are read without modification to obtain equivalent metrics when
   trimming is disabled. CheckM2 values come from `quality_report.tsv`.

   CleanGene automatically manages the CheckM2 DIAMOND database when
   `CHECKM2_MODE=required`, `CHECKM2_DB=""`, and
   `CHECKM2_AUTO_DOWNLOAD=true`. It checks the shared
   `<CleanGene>/databases/checkm2/` directory, reuses
   `CheckM2_database/uniref100.KO.1.dmnd` when present, and otherwise runs the
   installed CheckM2 database downloader once under a shared lock. Later runs
   reuse the same resolved database path. Set `CLEANGENE_DATABASE_ROOT` to place
   all managed databases on a dedicated shared filesystem; CheckM2 will use its
   `checkm2/` subdirectory. `CHECKM2_DATABASE_ROOT` can override only the
   CheckM2 parent, and `CHECKM2_DB` should be set only when intentionally using
   an exact custom `.dmnd` file.

   CleanGene resolves the CheckM2 executable at submission time, before the
   SLURM controller is started. Resolution honors `CHECKM2_EXECUTABLE`, then
   `PATH`, then an executable beside the active Python interpreter. The resolved
   path and version are written to run provenance and reused on resume. Run
   `cleangene doctor --config config/cleangene.arc.local.env` after editing
   local settings to verify the active checkout, executable, database roots, and
   scheduler access.

   Setting `CHECKM2_MODE=off` is allowed, but always produces `WARNING` with a
   Note that completeness and contamination were not evaluated. Do not turn it
   off merely to work around a database or environment problem; setup failures
   should be fixed and resumed.

   `CHECKM2_CPUS`, `CHECKM2_MEM`, and `CHECKM2_TIME` size the database setup
   job. Per-isolate CheckM2 prediction currently runs inside preprocess and uses
   `CPUS`; `CHECKM2_BATCH_SIZE` and `CHECKM2_MAX_INFLIGHT` are retained as
   compatibility settings and are not active scheduling controls.

   Global thresholds are configured with:

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
   warning-only. Numeric fail thresholds enable the normal PASS → WARNING →
   FAIL bands. Per-isolate manifest overrides use lowercase names, for example:

   ```text
   isolate_id  group_id  R1  R2  qc_max_contigs_pass  qc_min_coverage_fail
   ISO001      groupA    ... ... 250                  12
   ```

   Every global `QC_*` threshold shown above has a matching lowercase manifest
   column. The supported override columns are:

   ```text
   qc_max_contigs_pass
   qc_max_contigs_fail
   qc_min_n50_pass
   qc_min_n50_fail
   qc_min_coverage_pass
   qc_min_coverage_fail
   qc_min_read_length_pass
   qc_min_read_length_fail
   qc_min_mean_base_quality_pass
   qc_min_mean_base_quality_fail
   qc_min_completeness_pass
   qc_min_completeness_fail
   qc_max_checkm2_contamination_pass
   qc_max_checkm2_contamination_fail
   qc_max_kraken_contamination_fail
   ```

   Blank manifest cells inherit lower-priority values.

   `QC_PROFILE_FILE` accepts a tab-separated file with `scope_type`,
   `scope_value`, and any lowercase threshold columns:

   ```text
   scope_type  scope_value       qc_max_contigs_pass  qc_min_completeness_pass
   organism    Escherichia coli  250                  92
   group_id    outbreak_2026     200                  95
   ```

   `scope_type` must be `organism` or `group_id`. Resolution precedence is:

   ```text
   per-isolate manifest override
   > matching group_id profile
   > matching organism profile
   > global .env value
   > built-in default
   ```

   CleanGene rejects duplicate profiles, invalid/negative numbers, and invalid
   PASS/FAIL ordering. It copies the profile to `provenance/qc_profile.tsv` and
   writes every isolate's resolved thresholds to
   `provenance/qc_thresholds.tsv`, so resume does not depend on a changing
   external profile.

   Per-isolate and cohort `isolate_qc.tsv` outputs contain `PASS/FAIL`, `Notes`,
   `trimmed_read_length`, `mean_base_quality`, `sequencing_coverage`,
   `checkm2_completeness`, `checkm2_contamination`, and `qc_profile_source`.
   `PASS` means every evaluated criterion passed; `WARNING` means no hard
   failure but at least one warning or explicitly unavailable optional
   evaluation; `FAIL` means at least one hard failure. `Notes` contains all
   warnings and failures in deterministic semicolon-separated form. WARNING
   isolates remain in Panaroo and CleanGene read validation. Only FAIL isolates
   have `excluded=1` and are removed from downstream analysis.

> The script is written to work for many hundreds of isolates – the array ranges handle arbitrary sizes, and job dependencies ensure that stages start only after previous ones finish. Adjust defaults in `src/cleangene/defaults.py` if your resources differ.
