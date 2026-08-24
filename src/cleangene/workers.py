from __future__ import annotations
import csv, fcntl, gzip, hashlib, os, shutil, subprocess, sys, tempfile, time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from .config import assembler_mode, checkm2_mode, truthy
from .defaults import DEFAULTS
from .evidence import validate_isolate, validation_decision_logic_rows
from .fasta import assembly_metrics
from .manifest import groups, write_resolved
from .pangenome import normalize_panaroo, recover_sequences, select_rows, write_binary
from .plotting import plot_presence_absence
from .qc import QC_OUTPUT_FIELDS, classify_isolate_qc, isolate_thresholds, parse_checkm2_report, qc_value, read_metrics
from .slurm import array_task_count, assert_jobs_succeeded, available_slots, job_active, sbatch_cmd, submit_with_qos_retry, user_job_count, user_queue_snapshot
from .util import command_exists, load_json, read_tsv, run, safe_name, touch_done, write_tsv
from .ux import completed, log_line, waiting

STAGE_DESCRIPTIONS = (
    ("kraken_db_setup", "prepare the shared Kraken2 database when required"),
    ("preprocess", "check/trim adapters, run Kraken2 QC, optionally assemble with Shovill or direct SPAdes, and annotate with Prokka"),
    ("resolve_groups", "resolve organism groups and order the smallest groups first"),
    ("panaroo", "build and clean each group pangenome with Panaroo"),
    ("prepare_validation", "select genes and recover pangenome reference sequences"),
    ("validate", "map isolate reads with BWA and measure gene coverage, depth, and identity"),
    ("reduce", "apply read evidence and publish the final cleaned pangenome matrix"),
    ("plot", "render each group presence/absence summary"),
    ("summary", "compile cohort and group QC tables"),
)

ARRAY_STAGES = {"preprocess","panaroo","prepare_validation","validate","reduce","plot"}

def _controller_log(message: str, *, ok: bool = False) -> None:
    color=completed if ok else waiting
    print(color(log_line(message)),flush=True)

def context(run_dir: Path):
    cfg={**DEFAULTS,**load_json(run_dir/"provenance"/"resolved_config.json")}; rows=read_tsv(run_dir/"provenance"/"manifest.tsv"); return cfg, rows

def task_row(run_dir: Path, kind: str, index: int) -> dict[str,str]:
    rows=read_tsv(run_dir/"state"/f"{kind}_tasks.tsv")
    if index<0 or index>=len(rows): raise SystemExit(f"Task index {index} outside {kind} task list")
    return rows[index]

def manifest_row_for_task(task: dict[str,str], rows: list[dict[str,str]]) -> dict[str,str]:
    matches=[r for r in rows if r["isolate_id"]==task["isolate_id"]]
    if len(matches)!=1: raise SystemExit(f"Could not resolve manifest row for isolate {task['isolate_id']}")
    return {**matches[0], **task}

def user_excluded(row: dict[str,str]) -> bool:
    return truthy(row.get("user_excluded","false")) or truthy(row.get("exclude","false"))

def shlex_quote(x: str) -> str:
    import shlex
    return shlex.quote(x)

def parse_kraken_report(path: Path, expected: str) -> tuple[str,float,float]:
    top=("",-1.0); contamination=0.0; unclassified=0.0; expected_norm=" ".join(expected.casefold().split()); species=[]
    if not path.is_file(): return "",0.0,0.0
    for line in path.read_text(errors="replace").splitlines():
        f=line.split("\t")
        if len(f)<6: continue
        pct=float(f[0]); rank=f[3].strip(); name=" ".join(f[5].strip().split()); norm=name.lower()
        if rank=="S":
            species.append((name,pct))
            if pct>top[1]: top=(name,pct)
        if rank=="U": unclassified=max(unclassified,pct)
    reference=expected_norm or " ".join(top[0].casefold().split())
    contamination=sum(pct for name,pct in species if " ".join(name.casefold().split())!=reference)
    return top[0], contamination, max(0.0,100.0-unclassified)

def _run_checkm2(assembly: Path, out: Path, logs: Path, cfg: dict[str,str], isolate: str) -> tuple[float,float]:
    input_dir=out/"input"; result_dir=out/"results"; input_dir.mkdir(parents=True,exist_ok=True); result_dir.mkdir(parents=True,exist_ok=True)
    suffix=".fna.gz" if assembly.suffix==".gz" else ".fna"; link=input_dir/f"{safe_name(isolate)}{suffix}"
    _replace_symlink(link,assembly.resolve())
    run(["checkm2","predict","--threads",cfg.get("CPUS","4"),"--input",str(input_dir),"--output-directory",str(result_dir),"--database_path",cfg["CHECKM2_DB"],"--remove-intermediates","--force"],stdout=logs/"checkm2.stdout",stderr=logs/"checkm2.stderr")
    return parse_checkm2_report(result_dir/"quality_report.tsv")

def needs_kraken(rows: list[dict[str,str]], cfg: dict[str,str]) -> bool:
    mode=cfg.get("TAXONOMY_MODE","auto")
    return mode not in {"off","auto"} or any(r.get("grouping_source")=="kraken_pending" for r in rows)

def default_kraken2_db(run_dir: Path, cfg: dict[str,str]) -> Path:
    return Path(cfg.get("KRAKEN2_BUILD_DIR","") or run_dir.parent.parent/"databases"/f"kraken2_{cfg.get('KRAKEN2_DATABASE_SIZE','standard-8')}").expanduser().resolve()

def ensure_kraken2_db(run_dir: Path, cfg: dict[str,str], rows: list[dict[str,str]]) -> str:
    db=cfg.get("KRAKEN2_DB","").strip()
    if db and (Path(db)/"hash.k2d").is_file(): return db
    if not needs_kraken(rows,cfg):
        return ""
    if not truthy(cfg.get("KRAKEN2_AUTO_DOWNLOAD","true")):
        raise SystemExit("KRAKEN2_DB is required when taxonomy or Kraken-inferred grouping is enabled")
    db_path=default_kraken2_db(run_dir,cfg)
    if (db_path/"hash.k2d").is_file(): return str(db_path)
    raise SystemExit(f"Kraken2 database is not ready: {db_path}. The kraken_db_setup SLURM stage must complete before preprocess.")

def _directory_size(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())

def kraken_db_for_worker(db: str, cfg: dict[str,str]) -> tuple[str,bool]:
    mode=cfg.get("KRAKEN2_DB_ACCESS","auto").strip().lower()
    if mode not in {"auto","copy","mmap","direct"}:
        raise SystemExit("KRAKEN2_DB_ACCESS must be auto, copy, mmap, or direct")
    source=Path(db).resolve()
    if mode in {"mmap","direct"}: return str(source),mode=="mmap"
    cache_setting=cfg.get("KRAKEN2_NODE_CACHE_DIR","").strip()
    cache_root=Path(cache_setting).expanduser() if cache_setting else Path("/tmp")/f"cleangene-{os.getuid()}"/"kraken2"
    try:
        cache_root.mkdir(parents=True,exist_ok=True,mode=0o700)
        token=hashlib.sha256(f"{source}:{(source/'hash.k2d').stat().st_size}:{(source/'hash.k2d').stat().st_mtime_ns}".encode()).hexdigest()[:16]
        target=cache_root/f"{safe_name(source.name)}-{token}"
        lock_path=cache_root/f"{token}.lock"
        with lock_path.open("w") as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if not (target/".cleangene-ready").is_file():
                needed=_directory_size(source)+int(float(cfg.get("KRAKEN2_NODE_CACHE_MIN_FREE_GB","2"))*1024**3)
                if shutil.disk_usage(cache_root).free < needed: raise OSError(f"insufficient free space under {cache_root}")
                staging=Path(tempfile.mkdtemp(prefix=f"{token}.",dir=cache_root))
                try:
                    copied=staging/"db"; shutil.copytree(source,copied); (copied/".cleangene-ready").write_text(str(source)+"\n")
                    if target.exists(): shutil.rmtree(target)
                    os.replace(copied,target)
                finally:
                    shutil.rmtree(staging,ignore_errors=True)
        return str(target),True
    except OSError as error:
        if mode=="copy": raise SystemExit(f"Could not stage Kraken2 database in node-local cache {cache_root}: {error}")
        print(f"warning: Kraken2 node cache unavailable ({error}); using shared database with memory mapping",flush=True)
        return str(source),True

def _preprocess_scratch(cfg: dict[str,str], run_dir: Path, isolate: str) -> Path | None:
    if not truthy(cfg.get("PREPROCESS_USE_NODE_LOCAL_SCRATCH","true")): return None
    configured=cfg.get("PREPROCESS_SCRATCH_DIR","").strip()
    base=configured or os.environ.get("SLURM_TMPDIR","") or os.environ.get("TMPDIR","") or "/tmp"
    root=Path(base).expanduser()
    try:
        root.mkdir(parents=True,exist_ok=True)
        return Path(tempfile.mkdtemp(prefix=f"cleangene-{safe_name(run_dir.name)}-{safe_name(isolate)}-",dir=root))
    except OSError as error:
        print(f"warning: node-local preprocessing scratch unavailable ({error}); using run directory",flush=True)
        return None

def _sync_tree(source: Path, destination: Path) -> None:
    if source.is_dir(): shutil.copytree(source,destination,dirs_exist_ok=True)

def _shared_path(path: str, work_out: Path, shared_out: Path) -> str:
    try: return str(shared_out/Path(path).relative_to(work_out))
    except ValueError: return path

def kraken_db_setup(run_dir: Path, index: int | None = None) -> None:
    cfg, rows=context(run_dir)
    if not needs_kraken(rows,cfg): touch_done(run_dir/"state"/"kraken_db_setup.done.json",{"status":"not_required"}); return
    db=cfg.get("KRAKEN2_DB","").strip()
    db_path=Path(db).expanduser().resolve() if db else default_kraken2_db(run_dir,cfg)
    if not (db_path/"hash.k2d").is_file():
        db_path.parent.mkdir(parents=True,exist_ok=True)
        script=Path(__file__).resolve().parents[2]/"scripts"/"build_kraken2_database.sh"
        run([str(script),str(db_path),cfg.get("KRAKEN2_DB_CPUS",cfg.get("CPUS","4")),cfg.get("KRAKEN2_CLEAN_BUILD_FILES","true"),cfg.get("KRAKEN2_DATABASE_SIZE","standard-8")],stdout=run_dir/"logs"/"kraken2-db.stdout",stderr=run_dir/"logs"/"kraken2-db.stderr")
    cfg["KRAKEN2_DB"]=str(db_path)
    from .util import atomic_json
    atomic_json(run_dir/"provenance"/"resolved_config.json",cfg)
    touch_done(run_dir/"state"/"kraken_db_setup.done.json",{"KRAKEN2_DB":str(db_path)})

def _replace_symlink(link: Path, target: Path) -> None:
    link.parent.mkdir(parents=True,exist_ok=True)
    tmp=link.with_name(f".{link.name}.tmp-{os.getpid()}")
    if tmp.exists() or tmp.is_symlink(): tmp.unlink()
    os.symlink(str(target),tmp)
    os.replace(tmp,link)

def _link_manifest_fastqs(reads: Path, r1: str, r2: str) -> tuple[str,str]:
    link1=reads/"input_R1.fastq.gz"; link2=reads/"input_R2.fastq.gz"
    _replace_symlink(link1,Path(r1).expanduser().resolve())
    _replace_symlink(link2,Path(r2).expanduser().resolve())
    return str(link1),str(link2)

def _gzip_file(path: Path) -> Path:
    gz=path.with_name(path.name + ".gz")
    if gz.is_file(): return gz
    tmp=gz.with_name(f".{gz.name}.tmp-{os.getpid()}")
    with path.open("rb") as source, gzip.open(tmp,"wb",compresslevel=6) as target:
        shutil.copyfileobj(source,target)
    os.replace(tmp,gz)
    path.unlink()
    return gz

def _compress_assembly_outputs(assembly_dir: Path, assembly: Path, cfg: dict[str,str]) -> Path:
    mode=cfg.get("COMPRESS_ASSEMBLY_OUTPUTS","off").strip().lower()
    if mode not in {"off","intermediates","all"}:
        raise SystemExit("COMPRESS_ASSEMBLY_OUTPUTS must be off, intermediates, or all")
    if mode=="off" or not assembly_dir.is_dir(): return assembly
    for path in sorted(assembly_dir.rglob("*")):
        if not path.is_file() or path.suffix==".gz": continue
        if mode!="all" and path.absolute()==assembly.absolute(): continue
        if path.suffix.lower() in {".fa",".fasta",".gfa",".fastg"}: _gzip_file(path)
    gz_assembly=assembly.with_name(assembly.name + ".gz")
    return gz_assembly if mode=="all" and gz_assembly.is_file() else assembly

def _compress_annotation_outputs(annotation_dir: Path, gff: Path, cfg: dict[str,str]) -> None:
    mode=cfg.get("COMPRESS_ANNOTATION_OUTPUTS","off").strip().lower()
    if mode not in {"off","nonessential"}:
        raise SystemExit("COMPRESS_ANNOTATION_OUTPUTS must be off or nonessential")
    if mode=="off" or not annotation_dir.is_dir(): return
    for path in sorted(annotation_dir.iterdir()):
        if not path.is_file() or path.suffix==".gz" or path.resolve()==gff.resolve(): continue
        _gzip_file(path)

def prepare_read_inputs(row: dict[str,str], out: Path, logs: Path, cfg: dict[str,str]) -> tuple[str,str,str,int]:
    r1=row.get("R1","").strip(); r2=row.get("R2","").strip(); mode=cfg.get("READ_TRIMMING_MODE","auto").strip().lower()
    if truthy(cfg.get("SKIP_TRIM","false")) or assembler_mode(cfg) in {"spades","off"}: mode="off"
    method="manifest_fastq"; trimmed=0
    raw_bam=row.get("raw_bam","").strip()

    if raw_bam:
        if not Path(raw_bam).is_file():
            raise SystemExit(f"Input BAM not found for isolate {row.get('isolate_id','<unknown>')}: {raw_bam}")
        reads=out/"reads"; reads.mkdir(exist_ok=True)
        r1=str(reads/"raw_R1.fastq.gz"); r2=str(reads/"raw_R2.fastq.gz")
        collated=reads/"raw_name_collated.bam"
        other=reads/"raw_unpaired.fastq.gz"
        threads=str(max(1,int(cfg.get("CPUS","4"))-1))
        run(["samtools","collate","-@",threads,"-o",str(collated),raw_bam],stdout=logs/"samtools-collate.stdout",stderr=logs/"samtools-collate.stderr")
        run(["samtools","fastq","-@",threads,"-1",r1,"-2",r2,"-0",str(other),"-s",str(other),"-n",str(collated)],
            stdout=logs/"samtools-fastq.stdout",stderr=logs/"samtools-fastq.stderr")
        method="raw_bam_samtools_fastq"

    if not r1 or not r2:
        raise SystemExit(f"Missing R1/R2 read paths for isolate {row.get('isolate_id','<unknown>')}")

    missing=[p for p in (r1,r2) if not Path(p).is_file()]
    if missing:
        raise SystemExit(
            f"FASTQ input not found for isolate {row.get('isolate_id','<unknown>')}: "
            + ", ".join(missing)
        )
    if mode not in {"off","auto","always"}: raise SystemExit("READ_TRIMMING_MODE must be off, auto, or always")
    if mode=="off" and not raw_bam:
        r1,r2=_link_manifest_fastqs(out/"reads",r1,r2)
        method += "+symlinked"
    elif mode in {"auto","always"}:
        if command_exists("fastp"):
            reads=out/"reads"; reads.mkdir(exist_ok=True)
            tr1=str(reads/"trimmed_R1.fastq.gz"); tr2=str(reads/"trimmed_R2.fastq.gz")
            run(["fastp","--detect_adapter_for_pe","--in1",r1,"--in2",r2,"--out1",tr1,"--out2",tr2,"--thread",cfg.get("CPUS","4"),"--json",str(reads/"fastp.json"),"--html",str(reads/"fastp.html")],
                stdout=logs/"fastp.stdout",stderr=logs/"fastp.stderr")
            r1,r2=tr1,tr2
            method += "+fastp"
            trimmed=1
        elif mode=="always":
            raise SystemExit("READ_TRIMMING_MODE=always requires fastp")
        elif not raw_bam:
            r1,r2=_link_manifest_fastqs(out/"reads",r1,r2)
            method += "+symlinked"
    return r1,r2,method,trimmed

def preprocess(run_dir: Path, index: int) -> None:
    cfg, rows=context(run_dir); row=manifest_row_for_task(task_row(run_dir,"isolate",index),rows); iso=row["isolate_id"]; group=row["group_id"]; safe=safe_name(iso)
    root=run_dir/"results"/"groups"/safe_name(group); out=root/"01_isolates"/safe; done=run_dir/"state"/"preprocess"/f"{safe}.done.json"
    if done.is_file():
        status=load_json(done)
        if status.get("excluded") or status.get("external_pangenome") or (status.get("qc_status") and (out/"qc.tsv").is_file()) or (out/"annotation"/f"{safe}.gff").is_file(): return
    out.mkdir(parents=True,exist_ok=True); scratch=_preprocess_scratch(cfg,run_dir,iso); work_out=(scratch/"output") if scratch else out; work_out.mkdir(parents=True,exist_ok=True); logs=work_out/"logs"; logs.mkdir(exist_ok=True)
    fields=["isolate_id","group_id","excluded","reason","top_species","contamination_pct","R1","R2","raw_bam","read_preprocessing","adapter_trimmed","assembly","assembly_length","contigs","n50","l50","ambiguous_bases","gc_fraction","gff",*QC_OUTPUT_FIELDS]
    thresholds,profile_source=isolate_thresholds(run_dir,iso,cfg,row); mode=checkm2_mode(cfg)
    finished=False
    def finish(data: dict[str,object], payload: dict[str,object]) -> None:
        nonlocal finished
        if scratch: _sync_tree(work_out,out)
        for key in ("R1","R2","assembly","gff"):
            if data.get(key): data[key]=_shared_path(str(data[key]),work_out,out)
        if payload.get("gff"): payload["gff"]=_shared_path(str(payload["gff"]),work_out,out)
        write_tsv(out/"qc.tsv",fields,[data]); touch_done(done,payload); finished=True
    def assessment(*, expected: str = "", top: str = "", contamination: float | None = None,
                   reads: dict[str,float] | None = None, metrics: dict[str,object] | None = None,
                   completeness: float | None = None, checkm2_contamination: float | None = None,
                   internal: bool = False, external: bool = False, assembly: str = "",
                   gff_present: bool | None = None, warnings=(), errors=(), excluded: bool = False) -> dict[str,str]:
        metrics=metrics or {}; length=float(metrics["assembly_length"]) if metrics.get("assembly_length") not in {"",None} else None
        coverage=(reads["total_bases"]/length) if reads and length and length>0 else None
        return classify_isolate_qc(thresholds=thresholds,expected_organism=expected,top_species=top,
            kraken_contamination=contamination,read_length=None if not reads else reads["trimmed_read_length"],
            mean_quality=None if not reads else reads["mean_base_quality"],coverage=coverage,
            contigs=float(metrics["contigs"]) if metrics.get("contigs") not in {"",None} else None,
            n50=float(metrics["n50"]) if metrics.get("n50") not in {"",None} else None,
            completeness=completeness,checkm2_contamination=checkm2_contamination,checkm2_mode=mode,
            internal_pangenome=internal,external_pangenome=external,assembly_present=bool(assembly),
            gff_present=gff_present,user_exclusion=excluded,warnings=warnings,errors=errors)
    def qc_columns(result: dict[str,str], reads: dict[str,float] | None, metrics: dict[str,object] | None,
                   completeness: float | None, contamination: float | None) -> dict[str,object]:
        length=float(metrics["assembly_length"]) if metrics and metrics.get("assembly_length") not in {"",None} else None
        coverage=(reads["total_bases"]/length) if reads and length and length>0 else None
        return {**result,"trimmed_read_length":qc_value(None if not reads else reads["trimmed_read_length"]),
            "mean_base_quality":qc_value(None if not reads else reads["mean_base_quality"]),
            "sequencing_coverage":qc_value(coverage),"checkm2_completeness":qc_value(completeness),
            "checkm2_contamination":qc_value(contamination),"qc_profile_source":profile_source}
    try:
        if user_excluded(row):
            metrics={"assembly_length":"","contigs":"","n50":"","l50":"","ambiguous_bases":"","gc_fraction":""}
            result=assessment(excluded=True)
            finish({"isolate_id":iso,"group_id":group,"top_species":"","contamination_pct":"","R1":row.get("R1",""),"R2":row.get("R2",""),"raw_bam":row.get("raw_bam",""),"read_preprocessing":"skipped_user_excluded","adapter_trimmed":0,"assembly":"","gff":"",**metrics,**qc_columns(result,None,metrics,None,None)},{"excluded":True,"reason":result["reason"],"qc_status":result["PASS/FAIL"]}); return
        r1,r2,read_method,adapter_trimmed=prepare_read_inputs(row,work_out,logs,cfg)
        read_errors=[]
        try: read_qc=read_metrics(Path(r1),Path(r2),work_out/"reads"/"fastp.json")
        except (OSError,ValueError,SystemExit) as error: read_qc=None; read_errors.append(("read_metrics_failed",str(error)))
        expected=row.get("organism","").strip(); taxonomy=cfg.get("TAXONOMY_MODE","auto"); top=""; contam=None
        taxonomy_enabled=taxonomy not in {"off","auto"} or row.get("grouping_source")=="kraken_pending"
        if taxonomy_enabled:
            db=ensure_kraken2_db(run_dir,cfg,rows)
            if not db: raise SystemExit("KRAKEN2_DB is required when TAXONOMY_MODE is not off")
            worker_db,memory_map=kraken_db_for_worker(db,cfg); report=work_out/"kraken2.report.tsv"
            output=str(work_out/"kraken2.output.tsv") if truthy(cfg.get("KRAKEN2_KEEP_CLASSIFICATIONS","false")) else "/dev/null"
            command=["kraken2","--db",worker_db,"--threads",cfg.get("CPUS","4")]
            if memory_map: command.append("--memory-mapping")
            command += ["--paired","--report",str(report),"--output",output,r1,r2]
            run(command,stdout=logs/"kraken2.stdout",stderr=logs/"kraken2.stderr")
            top,contam,_=parse_kraken_report(report,expected)
        assembly=row.get("assembly","").strip()
        common={"isolate_id":iso,"group_id":group,"top_species":top,"contamination_pct":contam,"R1":r1,"R2":r2,"raw_bam":row.get("raw_bam",""),"read_preprocessing":read_method,"adapter_trimmed":adapter_trimmed}
        external_pangenome=bool(row.get("pangenome_dir","").strip())
        assembler=assembler_mode(cfg)
        internal_pangenome=not external_pangenome and assembler!="off"
        generated_assembly=False
        errors=list(read_errors); warnings=[]
        if not assembly and internal_pangenome:
            shov=work_out/"assembly"; shov.mkdir(exist_ok=True); assembly=str(shov/("contigs.fa" if assembler=="shovill" else "contigs.fasta"))
            generated_assembly=True
            if not Path(assembly).is_file():
                try:
                    if assembler=="shovill":
                        tmp=(scratch/"tmp"/"shovill") if scratch else out/"tmp"/"shovill"; tmp.mkdir(parents=True,exist_ok=True)
                        run(["shovill","--R1",r1,"--R2",r2,"--outdir",str(shov),"--tmpdir",str(tmp),"--cpus",cfg.get("CPUS","4"),"--force"],stdout=logs/"shovill.stdout",stderr=logs/"shovill.stderr")
                    else:
                        if any(shov.iterdir()): shutil.rmtree(shov)
                        run(["spades.py","--only-assembler","-1",r1,"-2",r2,"-o",str(shov),"-t",cfg.get("CPUS","4"),"-m",cfg.get("SPADES_MEMORY_GB","28")],stdout=logs/"spades.stdout",stderr=logs/"spades.stderr")
                except subprocess.CalledProcessError as error:
                    errors.append(("assembly_failed",f"{assembler} assembly failed with exit status {error.returncode}")); assembly=""
        if assembler=="off" and not assembly: warnings.append(("shovill_skipped","assembly and annotation were not evaluated because assembly mode is off"))
        empty_metrics={"assembly_length":"","contigs":"","n50":"","l50":"","ambiguous_bases":"","gc_fraction":""}
        if assembly and Path(assembly).is_file(): metrics=assembly_metrics(Path(assembly))
        else: metrics=empty_metrics; assembly=""
        completeness=check_contamination=None
        if mode=="required" and assembly:
            db=Path(cfg.get("CHECKM2_DB","")).expanduser()
            if not cfg.get("CHECKM2_DB","").strip() or not db.is_file(): errors.append(("checkm2_database_missing",f"CheckM2 database file was not found at '{cfg.get('CHECKM2_DB','')}'"))
            elif not command_exists("checkm2"): errors.append(("checkm2_unavailable","CheckM2 executable was not available"))
            else:
                try: completeness,check_contamination=_run_checkm2(Path(assembly),work_out/"checkm2",logs,cfg,iso)
                except (subprocess.CalledProcessError,ValueError,OSError) as error: errors.append(("checkm2_failed",f"CheckM2 evaluation failed: {error}"))
        pre=assessment(expected=expected,top=top,contamination=contam,reads=read_qc,metrics=metrics,
            completeness=completeness,checkm2_contamination=check_contamination,internal=internal_pangenome,
            external=external_pangenome,assembly=assembly,gff_present=None,warnings=warnings,errors=errors)
        if assembler=="off" and pre["PASS/FAIL"]!="FAIL": pre["reason"]="shovill_skipped"
        if not internal_pangenome:
            if generated_assembly and assembly: assembly=str(_compress_assembly_outputs(Path(assembly).parent,Path(assembly),cfg))
            data={**common,"assembly":assembly,"gff":"",**metrics,**qc_columns(pre,read_qc,metrics,completeness,check_contamination)}
            payload={"excluded":pre["PASS/FAIL"]=="FAIL","reason":pre["reason"],"qc_status":pre["PASS/FAIL"]}
            if external_pangenome: payload["external_pangenome"]=row["pangenome_dir"]
            finish(data,payload); return
        assembly_failures={"coverage_low","contigs_high","n50_low","completeness_low","checkm2_contamination_high","assembly_failed","assembly_missing","checkm2_database_missing","checkm2_unavailable","checkm2_failed"}
        if assembly_failures.intersection(pre["reason"].split(";")):
            skipped=assessment(expected=expected,top=top,contamination=contam,reads=read_qc,metrics=metrics,
                completeness=completeness,checkm2_contamination=check_contamination,internal=True,assembly=assembly,
                gff_present=False,warnings=warnings,errors=errors)
            if generated_assembly and assembly: assembly=str(_compress_assembly_outputs(Path(assembly).parent,Path(assembly),cfg))
            data={**common,"assembly":assembly,"gff":"",**metrics,**qc_columns(skipped,read_qc,metrics,completeness,check_contamination)}
            finish(data,{"excluded":True,"reason":skipped["reason"],"qc_status":"FAIL"}); return
        ann=work_out/"annotation"; gff=ann/f"{safe}.gff"; gff_present=True; prokka_errors=[]
        if not gff.is_file():
            if ann.exists(): shutil.rmtree(ann)
            try: run(["prokka","--outdir",str(ann),"--prefix",safe,"--locustag",safe,"--cpus",cfg.get("CPUS","4"),"--force",assembly],stdout=logs/"prokka.stdout",stderr=logs/"prokka.stderr")
            except subprocess.CalledProcessError as error: prokka_errors.append(("prokka_failed",f"Prokka failed with exit status {error.returncode}")); gff_present=None
        if gff_present is not None: gff_present=gff.is_file()
        final=assessment(expected=expected,top=top,contamination=contam,reads=read_qc,metrics=metrics,
            completeness=completeness,checkm2_contamination=check_contamination,internal=True,assembly=assembly,
            gff_present=gff_present,warnings=warnings,errors=[*errors,*prokka_errors])
        if gff.is_file(): _compress_annotation_outputs(ann,gff,cfg)
        if generated_assembly: assembly=str(_compress_assembly_outputs(Path(assembly).parent,Path(assembly),cfg))
        data={**common,"assembly":assembly,"gff":str(gff) if gff.is_file() else "",**metrics,**qc_columns(final,read_qc,metrics,completeness,check_contamination)}
        finish(data,{"excluded":final["PASS/FAIL"]=="FAIL","reason":final["reason"],"qc_status":final["PASS/FAIL"],"gff":str(gff) if gff.is_file() else ""})
    finally:
        if scratch:
            if not finished: _sync_tree(logs,out/"logs")
            shutil.rmtree(scratch,ignore_errors=True)

def retained_rows(run_dir: Path, group: str) -> list[dict[str,str]]:
    _, rows=context(run_dir); result=[]
    for row in rows:
        if row["group_id"]!=group: continue
        if user_excluded(row): continue
        safe=safe_name(row["isolate_id"]); qc=find_isolate_qc(run_dir,row)
        if not qc.is_file(): raise SystemExit(f"Missing isolate QC: {qc}")
        q=read_tsv(qc)[0]
        if q.get("PASS/FAIL")=="FAIL" or q["excluded"] in {"1","true","True"}: continue
        row=dict(row); row["assembly"]=q["assembly"]; row["gff"]=q["gff"]; row["R1"]=q.get("R1",row.get("R1","")); row["R2"]=q.get("R2",row.get("R2","")); result.append(row)
    return result

def find_isolate_qc(run_dir: Path, row: dict[str,str]) -> Path:
    safe=safe_name(row["isolate_id"])
    direct=run_dir/"results"/"groups"/safe_name(row["group_id"])/"01_isolates"/safe/"qc.tsv"
    if direct.is_file(): return direct
    hits=list((run_dir/"results"/"groups").glob(f"*/01_isolates/{safe}/qc.tsv"))
    return hits[0] if hits else direct

def _original_fastq_path(run_dir: Path, value: str) -> Path:
    """Resolve a manifest FASTQ without following a symlink at the final path."""
    path=Path(value).expanduser()
    if path.is_absolute(): return path
    candidates=[Path.cwd()/path]
    inputs=run_dir/"provenance"/"inputs.json"
    if inputs.is_file():
        manifest=load_json(inputs).get("manifest","")
        if manifest: candidates.append(Path(str(manifest)).expanduser().parent/path)
    return next((candidate.absolute() for candidate in candidates if candidate.is_file()),candidates[0].absolute())

def cleanup_trimmed_fastqs(run_dir: Path, dry_run: bool = False) -> dict[str,object]:
    """Replace retained fastp outputs with links to original manifest FASTQs."""
    _, rows=context(run_dir); results_root=(run_dir/"results").resolve(); report=[]; reclaimed=0
    for row in rows:
        iso=row["isolate_id"]
        if row.get("raw_bam","").strip():
            report.append([iso,"skipped_raw_bam","",0,"BAM-derived FASTQs must be retained for utilities"]); continue
        qc=find_isolate_qc(run_dir,row)
        if not qc.is_file():
            report.append([iso,"skipped_missing_qc","",0,str(qc)]); continue
        q=read_tsv(qc)[0]
        pairs=[]; problem=""
        for mate in ("R1","R2"):
            source=_original_fastq_path(run_dir,row.get(mate,"")); target=Path(q.get(mate,"")).expanduser()
            if not source.is_file(): problem=f"original {mate} not found: {source}"; break
            if not target.is_absolute(): target=(Path.cwd()/target).absolute()
            expected=f"trimmed_{mate}.fastq.gz"
            try: inside=target.parent.resolve().is_relative_to(results_root)
            except OSError: inside=False
            if target.name!=expected or not inside:
                problem=f"QC {mate} is not a CleanGene trimmed FASTQ: {target}"; break
            pairs.append((mate,source.absolute(),target))
        if problem:
            report.append([iso,"skipped",q.get("R1",""),0,problem]); continue
        sample_bytes=sum(target.lstat().st_size for _,_,target in pairs if target.exists() and not target.is_symlink())
        already=all(target.is_symlink() and target.resolve(strict=False)==source.resolve(strict=False) for _,source,target in pairs)
        if already:
            report.append([iso,"already_linked",str(pairs[0][2]),0,""]); continue
        if not dry_run:
            try:
                for _,source,target in pairs:
                    target.parent.mkdir(parents=True,exist_ok=True)
                    temporary=target.with_name(f".{target.name}.cleanup-{os.getpid()}")
                    if temporary.is_symlink() or temporary.exists(): temporary.unlink()
                    temporary.symlink_to(source)
                    os.replace(temporary,target)
            except OSError as error:
                if 'temporary' in locals() and (temporary.is_symlink() or temporary.exists()): temporary.unlink()
                report.append([iso,"error",str(pairs[0][2]),0,str(error)]); continue
        reclaimed += sample_bytes
        report.append([iso,"would_link" if dry_run else "linked",str(pairs[0][2]),sample_bytes,""])
    if not dry_run:
        write_tsv(run_dir/"results"/"cohort"/"fastq_cleanup.tsv",["isolate_id","status","trimmed_R1","bytes_reclaimed","detail"],report)
    counts={status:sum(1 for item in report if item[1]==status) for status in sorted({item[1] for item in report})}
    return {"isolates":len(report),"bytes_reclaimed":reclaimed,"counts":counts,"rows":report}

def compress_completed_outputs(run_dir: Path) -> dict[str,object]:
    """Compress safe run-local outputs, including preprocesses completed before resume."""
    cfg,rows=context(run_dir); assembly_mode=cfg.get("COMPRESS_ASSEMBLY_OUTPUTS","off").strip().lower(); annotation_mode=cfg.get("COMPRESS_ANNOTATION_OUTPUTS","off").strip().lower()
    if assembly_mode not in {"off","intermediates","all"}: raise SystemExit("COMPRESS_ASSEMBLY_OUTPUTS must be off, intermediates, or all")
    if annotation_mode not in {"off","nonessential"}: raise SystemExit("COMPRESS_ANNOTATION_OUTPUTS must be off or nonessential")
    report=[]; reclaimed=0
    for row in rows:
        qc=find_isolate_qc(run_dir,row)
        if not qc.is_file(): continue
        qc_rows=read_tsv(qc)
        if not qc_rows: continue
        q=qc_rows[0]; changed=False; iso=row["isolate_id"]
        assembly_dir=qc.parent/"assembly"
        if assembly_mode!="off" and assembly_dir.is_dir():
            recorded=Path(q.get("assembly","")).expanduser()
            recorded=recorded.absolute() if recorded else recorded
            for path in sorted(assembly_dir.rglob("*")):
                if not path.is_file() or path.is_symlink() or path.suffix==".gz": continue
                if assembly_mode!="all" and recorded and path.absolute()==recorded: continue
                if path.suffix.lower() not in {".fa",".fasta",".gfa",".fastg"}: continue
                before=path.stat().st_size; original=str(path); gz=_gzip_file(path); saved=max(0,before-gz.stat().st_size); reclaimed += saved
                report.append([iso,"assembly",original,str(gz),saved])
                if recorded and recorded==path.absolute(): q["assembly"]=str(gz); changed=True
        annotation_dir=qc.parent/"annotation"
        if annotation_mode=="nonessential" and annotation_dir.is_dir():
            gff=Path(q.get("gff","")).expanduser(); gff=gff.absolute() if gff else gff
            for path in sorted(annotation_dir.iterdir()):
                if not path.is_file() or path.is_symlink() or path.suffix==".gz" or path.suffix.lower()==".gff": continue
                if gff and path.absolute()==gff: continue
                before=path.stat().st_size; original=str(path); gz=_gzip_file(path); saved=max(0,before-gz.stat().st_size); reclaimed += saved
                report.append([iso,"annotation",original,str(gz),saved])
        if changed: write_tsv(qc,list(q.keys()),[q])
    write_tsv(run_dir/"results"/"cohort"/"storage_cleanup.tsv",["isolate_id","category","original_path","compressed_path","bytes_reclaimed"],report)
    counts={category:sum(1 for item in report if item[1]==category) for category in sorted({item[1] for item in report})}
    return {"files_compressed":len(report),"bytes_reclaimed":reclaimed,"counts":counts}

def group_size_class(n: int, cfg: dict[str,str]) -> str:
    if n <= int(cfg.get("PANAROO_SMALL_MAX_ISOLATES","499")): return "small"
    if n <= int(cfg.get("PANAROO_MEDIUM_MAX_ISOLATES","2000")): return "medium"
    return "large"

def build_organism_results_index(run_dir: Path) -> dict[str,int]:
    """Expose each complete isolate tree under its identified organism without copying data."""
    _,rows=context(run_dir); root=run_dir/"results"/"organisms"; root.mkdir(parents=True,exist_ok=True)
    records=[]; desired=set(); names={}
    for row in rows:
        qc=find_isolate_qc(run_dir,row)
        if not qc.is_file(): continue
        values=read_tsv(qc)
        if not values: continue
        q=values[0]; organism=(q.get("top_species","").strip() or row.get("organism","").strip()
            or (row.get("group_id","").strip() if not row.get("group_id","").startswith("__") else "") or "unclassified")
        organism_dir=root/safe_name(organism); link=organism_dir/safe_name(row["isolate_id"]); key=(organism_dir.name,link.name)
        previous=names.get(key)
        if previous and previous!=(organism,row["isolate_id"]):
            raise SystemExit(f"Organism results path collision: {previous} and {(organism,row['isolate_id'])}")
        names[key]=(organism,row["isolate_id"]); source=qc.parent.resolve(); organism_dir.mkdir(parents=True,exist_ok=True)
        if link.exists() and not link.is_symlink(): raise SystemExit(f"Refusing to replace non-symlink organism result path: {link}")
        relative=Path(os.path.relpath(source,link.parent)); _replace_symlink(link,relative); desired.add(link)
        records.append([organism,row["isolate_id"],str(link),str(source)])
    for path in root.glob("*/*"):
        if path.is_symlink() and path not in desired: path.unlink()
    for directory in root.iterdir():
        if directory.is_dir() and not directory.is_symlink() and not any(directory.iterdir()): directory.rmdir()
    write_tsv(run_dir/"results"/"cohort"/"organism_isolate_index.tsv",
        ["organism","isolate_id","symlink","target"],records)
    return {"organisms":len({row[0] for row in records}),"isolates":len(records)}

def resolve_groups(run_dir: Path, index: int | None = None) -> None:
    cfg, rows=context(run_dir)
    resolved=[]
    for row in rows:
        row=dict(row)
        if user_excluded(row):
            if row.get("grouping_source")=="kraken_pending":
                row["group_id"]="__user_excluded__"; row["grouping_source"]="user_excluded"
            resolved.append(row); continue
        if row.get("grouping_source")=="kraken_pending":
            qc=find_isolate_qc(run_dir,row)
            top=read_tsv(qc)[0].get("top_species","").strip() if qc.is_file() else ""
            if not top:
                raise SystemExit(f"Kraken2 could not infer organism for isolate {row['isolate_id']}")
            row["organism"]=top
            row["group_id"]=top
            row["grouping_source"]="kraken2_top_species"
        resolved.append(row)
    counts={g:sum(1 for r in resolved if r["group_id"]==g) for g in groups(resolved)}
    ordered=sorted(counts, key=lambda g:(counts[g],g))
    write_resolved(run_dir/"provenance"/"manifest.tsv",resolved)
    write_tsv(run_dir/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],([r["group_id"],r["isolate_id"]] for r in resolved))
    write_tsv(run_dir/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],([g,counts[g],group_size_class(counts[g],cfg)] for g in ordered))
    organism_index=build_organism_results_index(run_dir)
    touch_done(run_dir/"state"/"resolve_groups.done.json",{"groups":len(ordered),"order":"smallest_first","organism_index":organism_index})

def orchestrate_downstream(run_dir: Path, index: int | None = None) -> None:
    controller_downstream(run_dir)

def _stage_log_pattern(run_dir: Path, stage: str) -> Path:
    base=run_dir/"logs"/"slurm"
    return base/("preprocess" if stage=="preprocess" else "")/f"{stage}.%A_%a.log"

def _controller_cmd(run_dir: Path, cfg: dict[str,str], stage: str, array: str | None, cpus: str, mem: str, time_limit: str) -> list[str]:
    exe=f"{shlex_quote(sys.executable)} -m cleangene _worker"
    base=dict(account=cfg["SLURM_ACCOUNT"],partition=cfg["SLURM_PARTITION"])
    idx='${SLURM_ARRAY_TASK_ID}' if array else '0'
    wrap=f"{exe} --stage {stage} --run-dir {shlex_quote(str(run_dir))} --index {idx}"
    log=_stage_log_pattern(run_dir,stage); log.parent.mkdir(parents=True,exist_ok=True)
    return sbatch_cmd(name=f"cg-{stage}",wrap=wrap,cpus=cpus,mem=mem,time=time_limit,array=array,log=log,**base)

def _indices_spec(indices: list[int], max_parallel: str) -> str:
    if indices==list(range(indices[0],indices[-1]+1)):
        base=str(indices[0]) if len(indices)==1 else f"{indices[0]}-{indices[-1]}"
    else:
        base=",".join(map(str,indices))
    return f"{base}%{max_parallel}"

def _done(path: Path) -> bool:
    return path.is_file()

def incomplete_indices(run_dir: Path, stage: str) -> list[int]:
    if stage=="preprocess":
        rows=read_tsv(run_dir/"state"/"isolate_tasks.tsv")
        return [i for i,r in enumerate(rows) if not _done(run_dir/"state"/"preprocess"/f"{safe_name(r['isolate_id'])}.done.json")]
    rows=read_tsv(run_dir/"state"/"group_tasks.tsv")
    marker={"panaroo":"panaroo","prepare_validation":"prepare_validation","reduce":"reduce","plot":"plot"}[stage]
    return [i for i,r in enumerate(rows) if not _done(run_dir/"state"/marker/f"{safe_name(r['group_id'])}.done.json")]

def incomplete_validate_indices(run_dir: Path) -> list[int]:
    rows=read_tsv(run_dir/"state"/"isolate_tasks.tsv")
    return [i for i,r in enumerate(rows) if not _done(run_dir/"state"/"validate"/f"{safe_name(r['isolate_id'])}.done.json")]

def _float_or_none(value: str | None) -> float | None:
    try:
        return float(value) if value not in {"", "NA", None} else None
    except ValueError:
        return None

def _remove_done(path: Path) -> None:
    if path.is_file(): path.unlink()

def invalidate_legacy_identity_metrics(run_dir: Path, cfg: dict[str,str]) -> int:
    min_breadth=float(cfg["READ_VALIDATION_MIN_BREADTH"]); min_depth=float(cfg["READ_VALIDATION_MIN_MEAN_DEPTH"])
    invalidated=0
    for metrics in (run_dir/"results"/"groups").glob("*/03_read_validation/evidence/*/metrics.tsv"):
        rows=[]
        try:
            with metrics.open(newline="",errors="replace") as handle:
                rows=list(csv.DictReader(handle,delimiter="\t"))
        except OSError:
            continue
        stale=False
        for row in rows:
            mapped=_float_or_none(row.get("mapped_reads")) or 0
            breadth=_float_or_none(row.get("breadth")) or 0
            depth=_float_or_none(row.get("mean_depth")) or 0
            identity=_float_or_none(row.get("identity"))
            aligned=_float_or_none(row.get("aligned_positions"))
            method=row.get("identity_method","")
            if mapped>0 and breadth>=min_breadth and depth>=min_depth and identity==0 and (aligned in {0,None}) and not method:
                stale=True
                break
        if not stale: continue
        backup=metrics.with_name("metrics.pre_identity_fix.tsv")
        if not backup.is_file(): shutil.copy2(metrics,backup)
        group_safe=metrics.parents[3].name; iso_safe=metrics.parent.name
        _remove_done(run_dir/"state"/"validate"/f"{iso_safe}.done.json")
        _remove_done(run_dir/"state"/"reduce"/f"{group_safe}.done.json")
        _remove_done(run_dir/"state"/"plot"/f"{group_safe}.done.json")
        _remove_done(run_dir/"state"/"summary.done.json")
        invalidated += 1
    return invalidated

def invalidate_legacy_isolate_qc(run_dir: Path) -> int:
    invalidated=0; required=set(QC_OUTPUT_FIELDS)
    tasks=run_dir/"state"/"isolate_tasks.tsv"
    rows=read_tsv(tasks) if tasks.is_file() else []
    targets=[(run_dir/"state"/"preprocess"/f"{safe_name(row['isolate_id'])}.done.json",run_dir/"results"/"groups"/safe_name(row["group_id"])/"01_isolates"/safe_name(row["isolate_id"])/"qc.tsv") for row in rows]
    if not targets:
        targets=[(marker,hits[0] if hits else None) for marker in (run_dir/"state"/"preprocess").glob("*.done.json") for hits in [list((run_dir/"results"/"groups").glob(f"*/01_isolates/{marker.name[:-len('.done.json')]}/qc.tsv"))]]
    for marker,qc in targets:
        if not marker.is_file(): continue
        fields=set()
        if qc and qc.is_file():
            try:
                with qc.open(newline="",errors="replace") as handle: fields=set(next(csv.reader(handle,delimiter="\t"),[]))
            except OSError:
                fields=set()
        if required.issubset(fields): continue
        marker.unlink(); invalidated += 1
    if invalidated: _remove_done(run_dir/"state"/"summary.done.json")
    return invalidated

def resume_maintenance_signature(cfg: dict[str,str]) -> dict[str,object]:
    return {
        "legacy_identity_metrics":{
            "min_breadth":cfg.get("READ_VALIDATION_MIN_BREADTH",""),
            "min_depth":cfg.get("READ_VALIDATION_MIN_MEAN_DEPTH",""),
        },
        "isolate_qc_fields":list(QC_OUTPUT_FIELDS),
    }

def run_resume_maintenance(run_dir: Path, cfg: dict[str,str]) -> dict[str,int]:
    marker=run_dir/"state"/"resume_maintenance.done.json"; signature=resume_maintenance_signature(cfg)
    if marker.is_file():
        try:
            if load_json(marker).get("signature")==signature: return {"legacy_identity_metrics":0,"legacy_isolate_qc":0}
        except (OSError, ValueError):
            pass
    _controller_log("step=resume_maintenance | status=checking_legacy_validation_and_qc_markers")
    invalidated=invalidate_legacy_identity_metrics(run_dir,cfg)
    legacy_qc=invalidate_legacy_isolate_qc(run_dir)
    touch_done(marker,{"signature":signature,"legacy_identity_metrics":invalidated,"legacy_isolate_qc":legacy_qc})
    if invalidated: _controller_log(f"step=resume_maintenance | invalidated_legacy_identity_metrics={invalidated}",ok=True)
    if legacy_qc: _controller_log(f"step=resume_maintenance | invalidated_legacy_isolate_qc={legacy_qc}",ok=True)
    _controller_log("step=resume_maintenance | status=completed",ok=True)
    return {"legacy_identity_metrics":invalidated,"legacy_isolate_qc":legacy_qc}

def _wait_jobs(job_ids: list[str], cfg: dict[str,str], label: str, complete: str, details: str = "") -> None:
    poll=int(cfg["SLURM_POLL_SECONDS"])
    while True:
        active=job_active(job_ids)
        current=user_job_count()
        avail=available_slots(int(cfg["SLURM_USER_JOB_LIMIT"]),int(cfg["SLURM_JOB_HEADROOM"]),current)
        _controller_log(f"step={label} | user_jobs={current}/{cfg['SLURM_USER_JOB_LIMIT']} | available_slots={avail} | total_submitted={len(job_ids)} | total_completed=0 | current_step_completed={complete} | waiting_for_jobs")
        if not active:
            assert_jobs_succeeded(job_ids,details)
            return
        time.sleep(poll)

def _run_single_job(run_dir: Path, cfg: dict[str,str], stage: str, cpus: str, mem: str, time_limit: str, label: str) -> str:
    cmd=_controller_cmd(run_dir,cfg,stage,None,cpus,mem,time_limit)
    jid=submit_with_qos_retry(cmd,cfg,1,label)
    _wait_jobs([jid],cfg,label,"single job submitted",f"job_id={jid} stage={stage} index=0 log={_stage_log_pattern(run_dir,stage)}")
    return jid

def _index_done(run_dir: Path, stage: str, index: int) -> bool:
    kind="isolate" if stage in {"preprocess","validate"} else "group"
    row=task_row(run_dir,kind,index); name=row["isolate_id"] if kind=="isolate" else row["group_id"]
    marker={"validate":"validate"}.get(stage,stage)
    done_path=run_dir/"state"/marker/f"{safe_name(name)}.done.json"; done=_done(done_path)
    if stage=="reduce" and done:
        done=load_json(done_path).get("status")=="skipped" or _done(run_dir/"results"/"groups"/safe_name(name)/"cleaned_pangenome.tsv")
    return done

@dataclass
class _ActiveBatch:
    job_id: str
    stage: str
    indices: list[int]
    array: str
    seen: bool=False
    missing_polls: int=0

class _RollingScheduler:
    def __init__(self,run_dir: Path,cfg: dict[str,str]):
        self.run_dir=run_dir; self.cfg=cfg; self.active: dict[str,_ActiveBatch]={}; self.jobs=[]; self.submitted: dict[str,set[int]]={}; self.snapshot={"total":0,"jobs":{},"entries":[]}

    def _adopt_live_batches(self) -> None:
        discovered: dict[tuple[str,str],list[int]]={}
        run_arg=str(self.run_dir)
        for entry in self.snapshot.get("entries",[]):
            name=str(entry.get("name","")); command=str(entry.get("command",""))
            stage=name[3:] if name.startswith("cg-") else ""
            if stage not in ARRAY_STAGES or run_arg not in command: continue
            try: index=int(str(entry.get("task_id","")))
            except ValueError: continue
            discovered.setdefault((str(entry["job_id"]),stage),[]).append(index)
        for (jid,stage),indices in discovered.items():
            if jid not in self.active:
                unique=sorted(set(indices)); self.active[jid]=_ActiveBatch(jid,stage,unique,",".join(map(str,unique)),seen=True)
            self.submitted.setdefault(stage,set()).update(indices)

    def refresh(self) -> None:
        self.snapshot=user_queue_snapshot(); queued=self.snapshot["jobs"]
        self._adopt_live_batches()
        for jid,batch in list(self.active.items()):
            if jid in queued:
                batch.seen=True; batch.missing_polls=0; continue
            batch.missing_polls += 1
            if batch.missing_polls<2: continue
            details=f"job_id={jid} stage={batch.stage} index={batch.array} log={_stage_log_pattern(self.run_dir,batch.stage)}"
            assert_jobs_succeeded([jid],details)
            missing=[i for i in batch.indices if not _index_done(self.run_dir,batch.stage,i)]
            if missing: raise RuntimeError(f"SLURM job {jid} completed without {batch.stage} markers for indices: {','.join(map(str,missing[:10]))} | {details}")
            del self.active[jid]

    def active_indices(self,stage: str) -> set[int]:
        return {i for b in self.active.values() if b.stage==stage for i in b.indices}

    def stage_queue(self,stage: str) -> tuple[int,int]:
        running=pending=0; jobs=self.snapshot["jobs"]
        for jid,batch in self.active.items():
            if batch.stage!=stage: continue
            states=jobs.get(jid,{})
            running += states.get("RUNNING",0); pending += states.get("PENDING",0)+states.get("UNKNOWN",0)
        return running,pending

    def submit_ready(self,stage: str,indices: list[int],cpus: str,mem: str,time_limit: str,label: str,max_inflight: int) -> int:
        active_indices=self.active_indices(stage); ready=[i for i in indices if i not in active_indices and not _index_done(self.run_dir,stage,i)]
        submitted=0; max_batches=int(self.cfg["SLURM_MAX_OUTSTANDING_CHUNKS"]); chunk=int(self.cfg["SLURM_ARRAY_CHUNK_SIZE"])
        running,pending=self.stage_queue(stage); stage_queued=running+pending
        current=int(self.snapshot["total"]); avail=available_slots(int(self.cfg["SLURM_USER_JOB_LIMIT"]),int(self.cfg["SLURM_JOB_HEADROOM"]),current)
        batches=sum(1 for b in self.active.values() if b.stage==stage)
        while ready and avail>0 and stage_queued<max_inflight and batches<max_batches:
            n=min(chunk,avail,max_inflight-stage_queued,len(ready)); part=ready[:n]; ready=ready[n:]
            array=_indices_spec(part,self.cfg["SLURM_MAX_PARALLEL"]); cmd=_controller_cmd(self.run_dir,self.cfg,stage,array,cpus,mem,time_limit)
            jid=submit_with_qos_retry(cmd,self.cfg,array_task_count(array),label); batch=_ActiveBatch(jid,stage,part,array,seen=True)
            self.active[jid]=batch; self.jobs.append(jid); self.submitted.setdefault(stage,set()).update(part)
            local_jobs=self.snapshot["jobs"]; local_jobs[jid]={"PENDING":n}; self.snapshot["total"]=int(self.snapshot["total"])+n
            submitted += n; stage_queued += n; avail -= n; batches += 1
        return submitted

    def progress(self,stage: str,total_indices: list[int],label: str) -> None:
        complete=sum(_index_done(self.run_dir,stage,i) for i in total_indices); running,pending=self.stage_queue(stage)
        submitted=len({i for i in total_indices if _index_done(self.run_dir,stage,i)}|self.submitted.get(stage,set()))
        current=int(self.snapshot["total"]); avail=available_slots(int(self.cfg["SLURM_USER_JOB_LIMIT"]),int(self.cfg["SLURM_JOB_HEADROOM"]),current)
        total=len(total_indices); not_submitted=max(0,total-submitted)
        message=(
            f"step={label} | user_jobs={current}/{self.cfg['SLURM_USER_JOB_LIMIT']} | available_slots={avail} | "
            f"total_submitted={submitted} | total_completed={complete} | current_step_completed={complete}/{total} | "
            f"running={running} | slurm_pending={pending} | not_submitted_yet={not_submitted} | failed=0"
        )
        _controller_log(message)

    def wait_tick(self) -> None: time.sleep(int(self.cfg["SLURM_POLL_SECONDS"]))

def _stage_limit(cfg: dict[str,str],stage: str) -> int:
    if stage=="preprocess": return int(cfg["SLURM_PREPROCESS_MAX_INFLIGHT"])
    if stage=="validate": return int(cfg["SLURM_VALIDATION_MAX_INFLIGHT"])
    return int(cfg["SLURM_GROUP_MAX_INFLIGHT"])

def _run_index_stage(run_dir: Path, cfg: dict[str,str], stage: str, indices: list[int], cpus: str, mem: str, time_limit: str, label: str) -> list[str]:
    scheduler=_RollingScheduler(run_dir,cfg)
    while any(not _index_done(run_dir,stage,i) for i in indices):
        scheduler.refresh(); scheduler.submit_ready(stage,indices,cpus,mem,time_limit,label,_stage_limit(cfg,stage)); scheduler.progress(stage,indices,label)
        if any(not _index_done(run_dir,stage,i) for i in indices): scheduler.wait_tick()
    return scheduler.jobs

def _group_resources(cfg: dict[str,str],group_row: dict[str,str]) -> tuple[str,str,str]:
    prefix=f"PANAROO_{group_row.get('group_size_class','small').upper()}"
    return cfg.get(f"{prefix}_CPUS",cfg["PANAROO_CPUS"]),cfg.get(f"{prefix}_MEM",cfg["PANAROO_MEM"]),cfg.get(f"{prefix}_TIME",cfg["PANAROO_TIME"])

def _controller_pipeline(run_dir: Path,include_preprocess: bool) -> None:
    cfg,_=context(run_dir); group_rows=read_tsv(run_dir/"state"/"group_tasks.tsv"); isolate_rows=read_tsv(run_dir/"state"/"isolate_tasks.tsv"); scheduler=_RollingScheduler(run_dir,cfg)
    group_index={r["group_id"]:i for i,r in enumerate(group_rows)}; members={g:[] for g in group_index}
    for i,row in enumerate(isolate_rows): members.setdefault(row["group_id"],[]).append(i)
    rank={r["group_id"]:i for i,r in enumerate(group_rows)}; prep_order=sorted(range(len(isolate_rows)),key=lambda i:(rank.get(isolate_rows[i]["group_id"],len(rank)),i))
    all_prep=list(range(len(isolate_rows))); all_groups=list(range(len(group_rows))); all_validate=list(range(len(isolate_rows)))
    while True:
        scheduler.refresh()
        active=lambda stage:scheduler.active_indices(stage)
        group_prepared=lambda gi:all(_index_done(run_dir,"preprocess",i) for i in members.get(group_rows[gi]["group_id"],[]))
        panaroo_ready=[i for i in all_groups if group_prepared(i) and i not in active("panaroo") and not _index_done(run_dir,"panaroo",i)]
        for klass in ("small","medium","large"):
            ready=[i for i in panaroo_ready if group_rows[i].get("group_size_class")==klass]
            if ready:
                cpus,mem,limit=_group_resources(cfg,group_rows[ready[0]]); scheduler.submit_ready("panaroo",ready,cpus,mem,limit,"CleanGene panaroo",_stage_limit(cfg,"panaroo"))
        prepare_ready=[i for i in all_groups if _index_done(run_dir,"panaroo",i) and i not in active("prepare_validation") and not _index_done(run_dir,"prepare_validation",i)]
        for klass in ("small","medium","large"):
            ready=[i for i in prepare_ready if group_rows[i].get("group_size_class")==klass]
            if ready:
                cpus,mem,limit=_group_resources(cfg,group_rows[ready[0]]); scheduler.submit_ready("prepare_validation",ready,cpus,mem,limit,"CleanGene prepare_validation",_stage_limit(cfg,"prepare_validation"))
        validation_ready=[i for i,row in enumerate(isolate_rows) if _index_done(run_dir,"prepare_validation",group_index[row["group_id"]])]
        scheduler.submit_ready("validate",validation_ready,cfg["VALIDATION_CPUS"],cfg["VALIDATION_MEM"],cfg["VALIDATION_TIME"],"CleanGene validate",_stage_limit(cfg,"validate"))
        reduce_ready=[gi for gi,row in enumerate(group_rows) if _index_done(run_dir,"prepare_validation",gi) and all(_index_done(run_dir,"validate",i) for i in members[row["group_id"]])]
        scheduler.submit_ready("reduce",reduce_ready,cfg["SUMMARY_CPUS"],cfg["SUMMARY_MEM"],cfg["SUMMARY_TIME"],"CleanGene reduce",_stage_limit(cfg,"reduce"))
        plot_ready=[i for i in all_groups if _index_done(run_dir,"reduce",i)]
        scheduler.submit_ready("plot",plot_ready,cfg["PLOT_CPUS"],cfg["PLOT_MEM"],cfg["PLOT_TIME"],"CleanGene plot",_stage_limit(cfg,"plot"))
        if include_preprocess: scheduler.submit_ready("preprocess",prep_order,cfg["SLURM_CPUS"],cfg["SLURM_MEM"],cfg["SLURM_TIME"],"CleanGene preprocess",_stage_limit(cfg,"preprocess"))
        if include_preprocess and any(not _index_done(run_dir,"preprocess",i) for i in all_prep): scheduler.progress("preprocess",all_prep,"CleanGene preprocess")
        elif any(not _index_done(run_dir,"validate",i) for i in all_validate): scheduler.progress("validate",all_validate,"CleanGene validate")
        elif any(not _index_done(run_dir,"plot",i) for i in all_groups): scheduler.progress("plot",all_groups,"CleanGene group completion")
        pipeline_done=(not include_preprocess or all(_index_done(run_dir,"preprocess",i) for i in all_prep)) and all(_index_done(run_dir,"plot",i) for i in all_groups)
        if pipeline_done and not scheduler.active: break
        scheduler.wait_tick()
    if not _done(run_dir/"state"/"summary.done.json"): _run_single_job(run_dir,cfg,"summary",cfg["SUMMARY_CPUS"],cfg["SUMMARY_MEM"],cfg["SUMMARY_TIME"],"CleanGene summary")
    touch_done(run_dir/"state"/"orchestrate_downstream.done.json",{"groups":len(group_rows),"isolates":len(isolate_rows)})

def controller_downstream(run_dir: Path) -> None:
    cfg,_=context(run_dir); group_rows=read_tsv(run_dir/"state"/"group_tasks.tsv"); isolate_rows=read_tsv(run_dir/"state"/"isolate_tasks.tsv")
    for stage in ("panaroo","prepare_validation"):
        for klass in ("small","medium","large"):
            missing=set(incomplete_indices(run_dir,stage)); indices=[i for i,r in enumerate(group_rows) if r.get("group_size_class")==klass and i in missing]
            if indices:
                cpus,mem,limit=_group_resources(cfg,group_rows[indices[0]]); _run_index_stage(run_dir,cfg,stage,indices,cpus,mem,limit,f"CleanGene {stage}")
    val_indices=incomplete_validate_indices(run_dir)
    if val_indices: _run_index_stage(run_dir,cfg,"validate",val_indices,cfg["VALIDATION_CPUS"],cfg["VALIDATION_MEM"],cfg["VALIDATION_TIME"],"CleanGene validate")
    for stage,cpus,mem,limit in (("reduce",cfg["SUMMARY_CPUS"],cfg["SUMMARY_MEM"],cfg["SUMMARY_TIME"]),("plot",cfg["PLOT_CPUS"],cfg["PLOT_MEM"],cfg["PLOT_TIME"])):
        indices=incomplete_indices(run_dir,stage)
        if indices: _run_index_stage(run_dir,cfg,stage,indices,cpus,mem,limit,f"CleanGene {stage}")
    if not _done(run_dir/"state"/"summary.done.json"): _run_single_job(run_dir,cfg,"summary",cfg["SUMMARY_CPUS"],cfg["SUMMARY_MEM"],cfg["SUMMARY_TIME"],"CleanGene summary")
    touch_done(run_dir/"state"/"orchestrate_downstream.done.json",{"groups":len(group_rows),"isolates":len(isolate_rows)})

@contextmanager
def _controller_lock(run_dir: Path):
    path=run_dir/"state"/"controller.lock"; path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w") as handle:
        try: fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise SystemExit(f"Another CleanGene controller is already active for {run_dir}")
        handle.write(f"job_id={os.environ.get('SLURM_JOB_ID','unknown')} pid={os.getpid()}\n"); handle.flush(); yield

def slurm_controller(run_dir: Path, index: int | None = None) -> None:
    with _controller_lock(run_dir):
        cfg, rows=context(run_dir)
        _controller_log(f"controller_started | run_dir={run_dir} | isolates={len(rows)} | stages=" + " -> ".join(stage for stage,_ in STAGE_DESCRIPTIONS))
        run_resume_maintenance(run_dir,cfg)
        if needs_kraken(rows,cfg) and not cfg.get("KRAKEN2_DB","").strip() and not _done(run_dir/"state"/"kraken_db_setup.done.json"):
            _run_single_job(run_dir,cfg,"kraken_db_setup",cfg["KRAKEN2_DB_CPUS"],cfg["KRAKEN2_DB_MEM"],cfg["KRAKEN2_DB_TIME"],"CleanGene kraken_db_setup")
            cfg,rows=context(run_dir)
        unresolved=any(r.get("grouping_source")=="kraken_pending" for r in rows)
        if not unresolved and not _done(run_dir/"state"/"resolve_groups.done.json"):
            _run_single_job(run_dir,cfg,"resolve_groups",cfg["GROUP_ORCHESTRATOR_CPUS"],cfg["GROUP_ORCHESTRATOR_MEM"],cfg["GROUP_ORCHESTRATOR_TIME"],"CleanGene resolve_groups")
        if unresolved:
            prep_indices=incomplete_indices(run_dir,"preprocess")
            if prep_indices: _run_index_stage(run_dir,cfg,"preprocess",prep_indices,cfg["SLURM_CPUS"],cfg["SLURM_MEM"],cfg["SLURM_TIME"],"CleanGene preprocess")
            if not _done(run_dir/"state"/"resolve_groups.done.json"):
                _run_single_job(run_dir,cfg,"resolve_groups",cfg["GROUP_ORCHESTRATOR_CPUS"],cfg["GROUP_ORCHESTRATOR_MEM"],cfg["GROUP_ORCHESTRATOR_TIME"],"CleanGene resolve_groups")
            controller_downstream(run_dir)
        else:
            _controller_pipeline(run_dir,True)

def manifest_pangenome_dir(rows: list[dict[str,str]], group: str) -> Path | None:
    paths=sorted({r.get("pangenome_dir","").strip() for r in rows if r["group_id"]==group and r.get("pangenome_dir","").strip()})
    if not paths: return None
    if len(paths)>1: raise SystemExit(f"Group {group} has multiple pangenome_dir values: {', '.join(paths[:3])}")
    path=Path(paths[0]).expanduser().resolve()
    if not path.is_dir(): raise SystemExit(f"Prebuilt pangenome directory not found for group {group}: {path}")
    if not (path/"gene_presence_absence.csv").is_file(): raise SystemExit(f"Prebuilt pangenome is missing gene_presence_absence.csv: {path}")
    return path

def prepared_pangenome_dir(run_dir: Path, group: str, root: Path) -> Path:
    local=root/"02_pangenome"/"panaroo"
    if (local/"gene_presence_absence.csv").is_file(): return local
    done=run_dir/"state"/"panaroo"/f"{safe_name(group)}.done.json"
    if done.is_file():
        external=load_json(done).get("external_pangenome_dir","")
        if external: return Path(str(external))
    return local

def panaroo(run_dir: Path, index: int) -> None:
    cfg,rows=context(run_dir); task=task_row(run_dir,"group",index); group=task["group_id"]; root=run_dir/"results"/"groups"/safe_name(group); out=root/"02_pangenome"/"panaroo"; done=run_dir/"state"/"panaroo"/f"{safe_name(group)}.done.json"
    root.mkdir(parents=True,exist_ok=True); (root/"logs").mkdir(parents=True,exist_ok=True); out.mkdir(parents=True,exist_ok=True)
    if done.is_file():
        status=load_json(done)
        if status.get("status")=="skipped" or status.get("external_pangenome_dir") or (out/"gene_presence_absence.csv").is_file(): return
    retained=retained_rows(run_dir,group)
    if len(retained)<2: touch_done(done,{"status":"skipped","reason":"fewer_than_two_retained"}); return
    external=manifest_pangenome_dir(rows,group)
    logs=root/"logs"; logs.mkdir(parents=True,exist_ok=True)
    if external:
        isolates=[r["isolate_id"] for r in retained]; rows=normalize_panaroo(external/"gene_presence_absence.csv",isolates); calls=root/"02_pangenome"/"initial_calls"; calls.mkdir(parents=True,exist_ok=True); write_binary(calls/"gene_presence_absence.binary.tsv",rows,isolates)
        touch_done(done,{"n_isolates":len(isolates),"n_genes":len(rows),"external_pangenome_dir":str(external)})
        return
    gffs=[r["gff"] for r in retained if r.get("gff")]
    if len(gffs)<2: touch_done(done,{"status":"skipped","reason":"fewer_than_two_gffs"}); return
    threads=os.environ.get("SLURM_CPUS_PER_TASK",cfg.get("PANAROO_CPUS",cfg.get("CPUS","4")))
    run(["panaroo","-i",*gffs,"-o",str(out),"--clean-mode",cfg["PANAROO_CLEAN_MODE"],"-t",threads],stdout=logs/"panaroo.stdout",stderr=logs/"panaroo.stderr")
    isolates=[r["isolate_id"] for r in retained]; rows=normalize_panaroo(out/"gene_presence_absence.csv",isolates); calls=root/"02_pangenome"/"initial_calls"; calls.mkdir(parents=True,exist_ok=True); write_binary(calls/"gene_presence_absence.binary.tsv",rows,isolates)
    touch_done(done,{"n_isolates":len(isolates),"n_genes":len(rows)})

def prepare_validation(run_dir: Path, index: int) -> None:
    cfg,_=context(run_dir); group=task_row(run_dir,"group",index)["group_id"]; root=run_dir/"results"/"groups"/safe_name(group); done=run_dir/"state"/"prepare_validation"/f"{safe_name(group)}.done.json"
    root.mkdir(parents=True,exist_ok=True); (root/"logs").mkdir(parents=True,exist_ok=True)
    if done.is_file(): return
    retained=retained_rows(run_dir,group)
    if len(retained)<2: touch_done(done,{"status":"skipped"}); return
    isolates=[r["isolate_id"] for r in retained]; initial=root/"02_pangenome"/"initial_calls"/"gene_presence_absence.binary.tsv"; panaroo_dir=prepared_pangenome_dir(run_dir,group,root)
    if not initial.is_file(): touch_done(done,{"status":"skipped","reason":"missing_initial_pangenome"}); return
    rows=[]
    with initial.open(newline="") as h:
        for r in csv.DictReader(h,delimiter="\t"): rows.append({"Gene":r["Gene"],**{i:int(r[i]) for i in isolates}})
    selected=select_rows(rows,isolates,cfg["VALIDATION_SCOPE"],float(cfg["ACCESSORY_PREVALENCE_CUTOFF"])); out=root/"03_read_validation"; out.mkdir(parents=True,exist_ok=True)
    records,sources=recover_sequences(selected,panaroo_dir); from .fasta import write_fasta; write_fasta(out/"tested_gene_references.fasta",records); write_tsv(out/"tested_gene_key.tsv",["reference_id","Gene","sequence_source","source_locus","length"],sources)
    write_tsv(out/"selected_genes.tsv",["Gene"],([r["Gene"]] for r in selected))
    if records: run(["bwa","index",str(out/"tested_gene_references.fasta")],stdout=root/"logs"/"bwa-index.stdout",stderr=root/"logs"/"bwa-index.stderr")
    touch_done(done,{"n_selected":len(selected)})

def validate(run_dir: Path, index: int) -> None:
    cfg,_=context(run_dir); row=task_row(run_dir,"isolate",index); group=row["group_id"]; iso=row["isolate_id"]; root=run_dir/"results"/"groups"/safe_name(group); out=root/"03_read_validation"; ev=out/"evidence"/safe_name(iso); done=run_dir/"state"/"validate"/f"{safe_name(iso)}.done.json"
    if done.is_file():
        status=load_json(done)
        if status.get("status")=="skipped_excluded" or (ev/"metrics.tsv").is_file(): return
    retained={r["isolate_id"]:r for r in retained_rows(run_dir,group)}
    if iso not in retained: touch_done(done,{"status":"skipped_excluded"}); return
    key=out/"tested_gene_key.tsv"; ref=out/"tested_gene_references.fasta"
    if not key.is_file() or key.stat().st_size==0 or len(read_tsv(key))==0:
        ev.mkdir(parents=True,exist_ok=True); write_tsv(ev/"metrics.tsv",["reference_id","Gene","validated_call","validation_state","decision_reason","final_call_source","breadth","percent_coverage","mean_depth","identity","percent_identity","identity_method","identical_positions","aligned_positions","mapped_reads","mean_mapping_quality"],[]); touch_done(done,{"status":"no_selected_genes"}); return
    rr=retained[iso]; validate_isolate(ref,key,rr["R1"],rr["R2"],ev,int(cfg["VALIDATION_CPUS"]),float(cfg["READ_VALIDATION_MIN_BREADTH"]),float(cfg["READ_VALIDATION_MIN_MEAN_DEPTH"]),float(cfg["READ_VALIDATION_MIN_IDENTITY"]),int(cfg["READ_VALIDATION_MIN_MAPQ"]),int(cfg["BASEQUAL"])); touch_done(done)

def reduce_group(run_dir: Path, index: int) -> None:
    cfg,_=context(run_dir); group=task_row(run_dir,"group",index)["group_id"]; root=run_dir/"results"/"groups"/safe_name(group); done=run_dir/"state"/"reduce"/f"{safe_name(group)}.done.json"; cleaned=root/"cleaned_pangenome.tsv"
    if done.is_file() and cleaned.is_file(): return
    retained=retained_rows(run_dir,group)
    if len(retained)<2: touch_done(done,{"status":"skipped"}); return
    isolates=[r["isolate_id"] for r in retained]; initial_path=root/"02_pangenome"/"initial_calls"/"gene_presence_absence.binary.tsv"
    if not initial_path.is_file(): touch_done(done,{"status":"skipped","reason":"missing_initial_pangenome"}); return
    initial=read_tsv(initial_path); by_gene={r["Gene"]:{i:int(r[i]) for i in isolates} for r in initial}; out=root/"03_read_validation"; metrics=[]
    for iso in isolates:
        p=out/"evidence"/safe_name(iso)/"metrics.tsv"
        for r in read_tsv(p) if p.is_file() else []:
            r={**r,"isolate_id":iso,"initial_call":by_gene[r["Gene"]][iso]}
            r.setdefault("percent_coverage", str(float(r.get("breadth") or 0) * 100.0))
            if r.get("identity") not in {"", "NA", None}: r.setdefault("percent_identity", str(float(r.get("identity") or 0) * 100.0))
            else: r.setdefault("percent_identity","NA")
            r.setdefault("identity_method","legacy")
            r.setdefault("final_call_source","read_validation")
            r.setdefault("decision_reason","legacy metrics")
            metrics.append(r)
    validated={g:dict(v) for g,v in by_gene.items()}
    for r in metrics:
        if r.get("final_call_source")=="initial_call_unresolved":
            validated[r["Gene"]][r["isolate_id"]]=int(r["initial_call"])
        elif str(r.get("validated_call",""))!="":
            validated[r["Gene"]][r["isolate_id"]]=int(r["validated_call"])
    fields=["Gene",*isolates]; validated_matrix=out/"validated_gene_presence_absence.binary.tsv"; write_tsv(validated_matrix,fields,([g,*[validated[g][i] for i in isolates]] for g in by_gene)); shutil.copy2(validated_matrix,cleaned)
    metric_fields=["Gene","isolate_id","initial_call","validated_call","validation_state","decision_reason","final_call_source","breadth","percent_coverage","mean_depth","identity","percent_identity","identity_method","identical_positions","aligned_positions","mapped_reads","mean_mapping_quality"]
    write_tsv(out/"read_validation_metrics.tsv",metric_fields,metrics)
    evidence_rows=[]
    metric_index={(r["Gene"],r["isolate_id"]):r for r in metrics}
    for gene in by_gene:
        for iso in isolates:
            m=metric_index.get((gene,iso)); initial_call=by_gene[gene][iso]; final=validated[gene][iso]
            evidence_rows.append([group,gene,iso,initial_call,final,m["validation_state"] if m else "not_tested_carried_forward",m.get("final_call_source","not_tested") if m else "not_tested",m.get("breadth","") if m else "",m.get("percent_coverage","") if m else "",m.get("mean_depth","") if m else "",m.get("identity","") if m else "",m.get("percent_identity","") if m else "",m.get("identity_method","") if m else "",m.get("mapped_reads","") if m else ""])
    write_tsv(out/"gene_call_evidence.long.tsv",["group_id","Gene","isolate_id","initial_call","validated_call","validation_state","final_call_source","breadth","percent_coverage","mean_depth","identity","percent_identity","identity_method","mapped_reads"],evidence_rows)
    changes=[]
    for g in by_gene:
        a=[by_gene[g][i] for i in isolates]; b=[validated[g][i] for i in isolates]; changes.append([g,sum(a),sum(b),sum(x==0 and y==1 for x,y in zip(a,b)),sum(x==1 and y==0 for x,y in zip(a,b)),sum(x==y for x,y in zip(a,b))])
    write_tsv(out/"tested_genes.tsv",["Gene","n_initial_present","n_validated_present","n_added","n_removed","n_unchanged"],changes)
    touch_done(done,{"n_isolates":len(isolates),"n_genes":len(by_gene),"cleaned_pangenome":str(cleaned)})

def plot_group(run_dir: Path, index: int) -> None:
    cfg,_=context(run_dir); group=task_row(run_dir,"group",index)["group_id"]; root=run_dir/"results"/"groups"/safe_name(group); done=run_dir/"state"/"plot"/f"{safe_name(group)}.done.json"
    if done.is_file(): return
    matrix=root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv"; out=root/"04_summary"
    if not matrix.is_file():
        touch_done(done,{"status":"skipped","reason":"missing_validated_matrix"})
        return
    out.mkdir(parents=True,exist_ok=True)
    plot_presence_absence(matrix,out,group,int(cfg["PLOT_MAX_CLUSTER_ISOLATES"]))
    touch_done(done,{"matrix":str(matrix),"outdir":str(out)})

def summarize(run_dir: Path) -> None:
    cfg,rows=context(run_dir); storage=compress_completed_outputs(run_dir); iso_rows=[]; group_rows=[]
    for group in groups(rows):
        root=run_dir/"results"/"groups"/safe_name(group); retained=retained_rows(run_dir,group); val=root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv"
        n_genes=max(0,len(read_tsv(val))) if val.is_file() else 0; group_rows.append([group,len([r for r in rows if r["group_id"]==group]),len(retained),n_genes])
        for row in rows:
            if row["group_id"]!=group: continue
            q=read_tsv(find_isolate_qc(run_dir,row))[0]; q={**q,"group_id":group}
            if user_excluded(row): q["excluded"]="1"; q["reason"]="user_excluded"
            iso_rows.append(q)
    cohort=run_dir/"results"/"cohort"; write_tsv(cohort/"isolate_qc.tsv",["isolate_id","group_id","excluded","reason","top_species","contamination_pct","R1","R2","raw_bam","read_preprocessing","adapter_trimmed","assembly","assembly_length","contigs","n50","l50","ambiguous_bases","gc_fraction","gff",*QC_OUTPUT_FIELDS],iso_rows); write_tsv(cohort/"group_summary.tsv",["group_id","input_isolates","retained_isolates","validated_gene_clusters"],group_rows); write_tsv(cohort/"validation_decision_logic.tsv",["state","criteria","final_call_behavior","biological_interpretation"],validation_decision_logic_rows(cfg["READ_VALIDATION_MIN_BREADTH"],cfg["READ_VALIDATION_MIN_MEAN_DEPTH"],cfg["READ_VALIDATION_MIN_IDENTITY"])); payload={"groups":len(group_rows),"storage_cleanup":storage,"organism_index":build_organism_results_index(run_dir)}
    if truthy(cfg.get("CLEANUP_TRIMMED_FASTQ","false")):
        cleanup=cleanup_trimmed_fastqs(run_dir); payload["fastq_cleanup"]={key:value for key,value in cleanup.items() if key!="rows"}
    touch_done(run_dir/"state"/"summary.done.json",payload)

def dispatch(stage: str, run_dir: Path, index: int | None) -> None:
    if stage=="slurm_controller": slurm_controller(run_dir,index)
    elif stage=="kraken_db_setup": kraken_db_setup(run_dir,index)
    elif stage=="preprocess": preprocess(run_dir,int(index))
    elif stage=="resolve_groups": resolve_groups(run_dir,index)
    elif stage=="orchestrate_downstream": orchestrate_downstream(run_dir,index)
    elif stage=="panaroo": panaroo(run_dir,int(index))
    elif stage=="prepare_validation": prepare_validation(run_dir,int(index))
    elif stage=="validate": validate(run_dir,int(index))
    elif stage=="reduce": reduce_group(run_dir,int(index))
    elif stage=="plot": plot_group(run_dir,int(index))
    elif stage=="summary": summarize(run_dir)
    else: raise SystemExit(f"Unknown worker stage: {stage}")
