
## Quick Start on a SLURM cluster

1. **Install Conda/Rocks**   
   ```bash
   mamba env create -f environment.yml
   conda activate cleangene
   pip install -e .
   ```

2. **Prepare configuration**  
   Copy the example and fill your paths.

   ```bash
   cp config/cleangene.example.env config/cleangene.env
   vi config/cleangene.env
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
   Job status with `squeue -j <JOB_ID>`. Logs are in `<analysis_root>/runs/<run-id>/logs/slurm/*.log`.

7. **Large scale (≥30 k isolates)**  
   Use `config/cleangene.arc.env` on ARC. It leaves `SLURM_PARTITION` empty, submits organisms smallest first, and uses Panaroo resource tiers for small (<500), medium (500-2000), and large (>2000) organism groups.

8. **Resume an unfinished run**

   ```bash
   cleangene run \
       --analysis-root ~/cleangene-output \
       --config config/cleangene.env \
       --resume <run-id>
   ```

9. **Local debugging** (small tests)

   ```bash
   cleangene run \
       --manifest input/manifest.tsv \
       --analysis-root ./local_run \
       --profile local
   ```

> The script is written to work for many hundreds of isolates – the array ranges handle arbitrary sizes, and job dependencies ensure that stages start only after previous ones finish. Adjust defaults in `src/cleangene/defaults.py` if your resources differ.
