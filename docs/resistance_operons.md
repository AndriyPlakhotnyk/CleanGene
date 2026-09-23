# Resistance operon analysis

`--resistance-operon` adds a stage after cleaned pangenome reduction and plots,
before final cohort archiving. The requested spelling `-ressitanec-operon` is an
alias. It works with `run`, `run --resume`, and `resume`, locally and through Slurm.
It is disabled by default. Existing assembly, annotation, pangenome and read
validation behavior is retained.

## Submit the E. faecium cohort on ARC

Activate the updated CleanGene environment. `environment.yml` now includes
`ncbi-amrfinderplus`, `edlib`, and the existing MAFFT/mapping/plotting tools.
The analysis needs an installed AMRFinderPlus database; it does not let hundreds
of array tasks race to download or update one.

`AMRFINDER_DB` must be set to an existing, versioned shared database directory.
CleanGene verifies that exact directory with `amrfinder --database_version
--database <path>` during preflight, before submitting any analysis arrays. A
missing or invalid path stops the controller with AMRFinder's diagnostic output;
it is not silently replaced with a node-local default database.

For a fresh ARC checkout, or to fast-forward an existing clean checkout to the
latest `main`, use the bootstrap wrapper. It refuses to overwrite local Git
changes, installs both Conda environments, updates AMRFinderPlus in shared
storage, and records the resulting versioned database path:

```bash
module load git
module load miniforge3       # use the Miniforge module name provided by ARC
bash scripts/bootstrap_arc_latest.sh \
  /project/def-YOURLAB/cleangene \
  /project/def-YOURLAB/conda-envs \
  /project/def-YOURLAB/databases/amrfinderplus \
  main
```

The wrapper is safe to rerun after later updates. If ARC uses another module
name, load that module so `mamba` is on `PATH` first. See
[scripts/bootstrap_arc_latest.sh](../scripts/bootstrap_arc_latest.sh).

```bash
conda activate cleangene
# Update the environment first if AMRFinderPlus/edlib are missing:
# bash scripts/install_or_update.sh

# One-time database setup; use a writable project/shared storage path.
amrfinder_update --database /path/to/shared/databases/amrfinderplus

cp config/cleangene.arc.env config/cleangene.resistance.arc.local.env
```

Edit that local config to set your real account, partition and **versioned**
database directory, for example:

```ini
SLURM_ACCOUNT="YOUR_ARC_ACCOUNT"
SLURM_PARTITION="YOUR_ARC_PARTITION"
AMRFINDER_DB="/path/to/shared/databases/amrfinderplus/YYYY-MM-DD.N"
AMRFINDER_ORGANISM="Enterococcus_faecium"
```

Do not use a moving `latest` database link for a long production cohort. The tool
and database versions and actual commands are retained with the results. The
manifest should contain `isolate_id`, `group_id`, `R1`, `R2`, with
`Enterococcus faecium` as the group for this cohort and absolute input paths.
Taxonomy follows the normal ARC configuration: `TAXONOMY_MODE=auto` skips
Kraken2 for an explicitly supplied group. To also classify these known-species
inputs and assess contamination, set `TAXONOMY_MODE="contamination"` and configure
`KRAKEN2_DB` as for an ordinary run; add an `organism` manifest column with the
expected species. The resistance flag does not change taxonomic QC behavior.

```bash
# Preview the Slurm controller submission:
bash scripts/submit_resistance_arc.sh input/efaecium.manifest.tsv \
  config/cleangene.resistance.arc.local.env --run-id efaecium_resistance --dry-run

# Submit from the login node (the launcher itself calls sbatch):
bash scripts/submit_resistance_arc.sh input/efaecium.manifest.tsv \
  config/cleangene.resistance.arc.local.env --run-id efaecium_resistance
```

The launcher executes this full pipeline command:

```bash
python -m cleangene run \
  --manifest input/efaecium.manifest.tsv \
  --config config/cleangene.resistance.arc.local.env \
  --analysis-root "$PWD" \
  --profile slurm \
  --assembler shovill \
  --skip-downsampling \
  --ignore-checkm2 \
  --resistance-operon \
  --run-id efaecium_resistance
```

Shovill, annotation, Panaroo, read validation, arbitration, QC and figures still
run. Raw-read depth reduction and CheckM2 setup/prediction are disabled. Read
trimming retains its normal configuration. QC exclusions still apply. At around
1,000 isolates, the ARC template selects the medium Panaroo tier (32 CPUs,
160 GB, 72 hours), which fits ordinary ARC CPU nodes with headroom. If
`arc.hardware` and `sinfo` show that your account can use a big-memory
partition, you may raise `PANAROO_MEDIUM_MEM` in the local config; otherwise a
512 GB request is rejected before the array is submitted.
Resistance detection uses bounded arrays (100 concurrent tasks, 4 CPUs/16 GB
per task by default) and an 8-CPU/64-GB merge job. Array chunks are bounded by the
existing scheduler. No claim of a 1,000-isolate runtime benchmark is implied.

To add the stage to a completed run or resume a failed submission:

```bash
python -m cleangene resume --run-dir "$PWD/runs/efaecium_resistance" \
  --config config/cleangene.resistance.arc.local.env \
  --ignore-checkm2 --skip-downsampling --resistance-operon
```

## What is being called an operon

AMRFinderPlus detects resistance determinants; it does not determine
co-transcription or general operon boundaries. The outputs therefore describe
**candidate resistance operons or single-determinant genomic contexts**.
Combined nucleotide/protein/GFF searches are used for available Prokka proteins,
including compressed `.faa.gz`. Otherwise nucleotide search is recorded in the
command provenance. E. faecium-specific mutation screening is enabled for that
organism. Non-AMR plus hits are retained in raw output but excluded from these
resistance-locus summaries. Point-mutation hits are grouped by the underlying gene
(e.g. gyrA), while their original mutation symbols remain in `reported_amr_symbols`;
mutation alleles do not create spurious gene-presence categories.

Built-in `vanA` and `vanB` gene sets define candidate spans:

- vanA: vanR-A, vanS-A, vanH-A, vanA, vanX-A, vanY-A, vanZ-A.
- vanB: vanR-B, vanS-B, vanY-B, vanW-B, vanH-B, vanB, vanX-B.

These are explicit analysis definitions, not a claim that every natural locus
contains every listed gene or that every gene is in one transcriptional unit.
Known AMRFinder symbols and exact annotation gene names extend a span on the
**same contig**, subject to a maximum 10-kb gap between members. Other
determinants default to the AMR gene span. Every intervening CDS is included,
plus exactly one nearest nonoverlapping CDS on each side. The region is oriented
to the defining determinant; upstream/downstream roles follow that orientation.
Repeated nonoverlapping copies of the same determinant start separate loci.
Fragmented contigs are never joined speculatively. Missing expected members,
missing flanks and partial AMR hits remain explicit flags.

To change or add definitions, set `RESISTANCE_OPERON_DEFINITIONS` to an absolute
JSON path mapping type names to exact AMRFinder/annotation gene symbols. Definitions
with the same name override the built-ins, and a symbol cannot belong to two
sets. The resolved definitions are saved in provenance. For example, a reviewed
custom definition file could contain:

```json
{"my_reviewed_type": ["exact_gene_symbol_1", "exact_gene_symbol_2"]}
```

Review definitions and fragmentation flags before interpreting a locus as a
complete operon. A missing assembly segment is not evidence of gene deletion.

## Categories, bins and variants

1. **Gene-content category:** binary presence/absence of all CDS families inside
   the extracted region, including both flanking genes. AMR identities use
   AMRFinder symbols; other genes use per-isolate Panaroo family assignments.
   Genes unavailable in the Panaroo mapping are clustered at 95% global nucleotide
   similarity. Hypothetical proteins are never merged merely because they share
   the same product label. Missing flanks use explicit `UNKNOWN_*` markers.
2. **95% similarity bin:** within each gene-content category, compare the complete
   oriented nucleotide sequences, including intergenic DNA and flanks. Similarity
   is `1 - global edit distance / max(sequence lengths)`. SNPs and inserted/deleted
   bases count against similarity; ambiguous bases do not count as matches.
   The threshold is inclusive. Abundance-ordered greedy clustering assigns each
   member directly to a fixed representative; length and stable IDs break ties.
   This avoids single-linkage chains. Two members in the same bin need not be
   95% similar to one another; both meet the threshold against the representative.
3. **Exact variants:** retain each distinct assembly sequence within its category,
   even if it belongs to the same 95% bin. Count distinct isolates, retaining copy
   counts separately. Rank variants per resistance type, then display the top 12.
   All unique sequences, including rare variants beyond the top 12, are aligned
   with MAFFT FFT-NS-1 and receive SNP/indel event tables.

Each type uses its most prevalent exact variant as the alignment reference, with
stable ID tie-breaking. This is not an ancestral or susceptible reference. Events
use 1-based MSA columns and oriented locus coordinates; insertions use the
preceding reference base (0 before the first base). Indels below 50 bp are marked
small and those at least 50 bp large. Uncertain bases are distinguished from SNPs.
Gene copy counts, order, strand and length differences are also reported. An MSA
can display large/repeat-associated differences but does not establish an
inversion/transposition mechanism or a trustworthy breakpoint in repetitive DNA.

Categories and bins can coexist in one isolate, so their prevalence percentages
need not sum to 100%. Denominators include all retained, successfully evaluated
isolates in the organism group, including zero-hit isolates. Excluded isolates
are `NA`. No-hit means **not detected in this assembly**, not proven absence.
Partial and read-warning loci remain visible rather than silently disappearing.

## Read support and replicon origin

Variant sequences are assembly-derived. When reads are available, the stage maps
them to the **whole isolate assembly**, so homologous regions compete for mapping.
Default support requires MAPQ >=20, base quality >=30, depth >=5, and >=90%
assembly-allele agreement at a position. A locus is `assembly_and_reads` when
>=95% of positions pass and there are no sufficiently covered discordant
positions. Insertions/deletions in the pileup count against support. Low-MAPQ
repeats, low coverage, and allele disagreement produce `assembly_read_warning`;
no available reads or disabled read checks produce `assembly_only`.

Compressed per-base evidence is retained. Mapping BAMs are discarded by default
after evidence extraction to limit storage at cohort scale; set
`RESISTANCE_KEEP_BAM=true` to retain alignments (the normal final stage archives
them as verified CRAM files). Read-supported locus status does not
mean every nucleotide is resolved, nor does short-read mapping establish complete
molecular phasing, repeat copy number or expression. Use the per-base evidence
and contig context for critical variant review. There is no automatic minority
haplotype reconstruction or antibiotic susceptibility prediction.

Every locus includes `origin` and `origin_evidence`. The default is `unknown`.
Neither AMR gene identity nor contig size is used to guess chromosome/plasmid.
Supply reviewed contig-level evidence through `RESISTANCE_CONTIG_ORIGINS`, a TSV:

```text
isolate_id	contig	origin	evidence
sample1	contig00001	chromosome	closed reference-matched chromosome
sample1	contig00003	plasmid	MOB-suite prediction; review required
```

Allowed origins are `chromosome`, `plasmid`, `unknown`; evidence is mandatory and
preserved verbatim. These are supplied assignments, not classifications performed
by CleanGene. A partial contig without reliable origin evidence stays unknown.

## Outputs

All outputs are under `runs/<run>/results/resistance_analysis/`:

| Folder/file | Contents |
|---|---|
| `figures/operon_prevalence.*` | Cohort prevalence of detected resistance types |
| `figures/<organism>/<type>/prevalence.*` | Gene-content category and top exact-variant prevalence |
| `figures/<organism>/<type>/alignment_overview.*` | Top 12 variants with SNP/insertion/deletion/uncertainty colours |
| `figures/<organism>/<type>/alignment_sites.pdf` | Paginated base-resolution differences; individual SNPs remain visible |
| `figures/<organism>/<type>/gene_structure.*` | Gene order, orientation and spacing including flanks |
| `tables/loci.tsv` | Per-copy bins, coordinates, origin, completeness and read support |
| `tables/categories.tsv`, `similarity_bins.tsv`, `variants.tsv` | Three levels of grouping and distinct-isolate counts |
| `tables/variant_events.tsv`, `gene_structure.tsv` | All exact variants, including rare variants beyond the plotted 12 |
| `tables/locus_genes.tsv` | Per-gene coordinates, family, annotation and flank role |
| `tables/<organism>/<type>/gene_presence_absence.tsv` | Per-locus gene content including flanks and expected members |
| `tables/isolate_presence_absence.tsv`, `sample_status.tsv` | Detection status and denominator audit |
| `alignments/<organism>/<type>/` | All unique raw and aligned nucleotide FASTAs plus MAFFT log |
| `isolates/<sample>/` | Raw AMRFinderPlus output, commands, logs, locus data, compressed base support |
| `provenance/` | Settings, definitions and AMRFinder software/database version |

Resume checks isolate input metadata and analysis settings, reusing unchanged
scans. Missing required artifacts cause rescan; changed settings or input metadata
invalidate scans. Reports are regenerated when scan outputs, Panaroo assignments,
origin evidence or settings change, or when an output is missing. Generated
`tables`, `figures`, and `alignments` directories are replaced when regeneration
is needed; keep manual edits elsewhere.

## Validation and limits

Run the regression suite in the CleanGene environment:

```bash
bash tests/run_tests.sh
```

A separate **real-tool**, two-isolate execution test exercises Shovill, Prokka,
Panaroo, read validation, AMRFinderPlus, MAFFT and figures, with downsampling and
CheckM2 off. It checks an expected SNP and 3-bp indel, two exact variants in one
95% bin, read support and resume reuse:

```bash
PYTHONPATH=src python scripts/test_resistance_e2e.py \
  --work-dir test_data/resistance-e2e \
  --amrfinder-db /path/to/versioned/amrfinder/database \
  --genome /path/to/bacterial-reference.fasta
```

The fixture uses a 200-kb bacterial backbone with a public AMRFinder reference in
an in-silico sequence and generated reads. It tests software execution, not
accuracy on the forthcoming clinical E. faecium dataset. The test runs with taxonomic
QC disabled because its sequence is a computational fixture. Production retains
the chosen taxonomy configuration. This fixture does not validate Kraken2 database
classification accuracy. Slurm submission/resource/index behavior is
tested without a live ARC allocation; a full 1,000-isolate production run has not
been performed here.

AMRFinder interface and interpretation references:
[NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/),
[running AMRFinderPlus](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus),
[interpreting results](https://github.com/ncbi/amr/wiki/Interpreting-results).
