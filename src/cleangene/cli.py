from __future__ import annotations
import argparse, json, os, shutil, subprocess, sys, tempfile, time
import re
from datetime import datetime
from pathlib import Path
from .config import assembler_mode, checkm2_mode, read_env, truthy
from .checkm2 import CheckM2DbError, CheckM2DbNotReady, bundled_test_genome, checkm2_database_root, checkm2_named_input_link, checkm2_predict_capabilities, checkm2_predict_command, checkm2_testrun_command, parse_checkm2_quality_report, resolve_checkm2_db
from .completion import reconcile_preprocess_outputs
from .defaults import DEFAULTS, SCIENTIFIC_DEFAULTS
from .manifest import groups, load_manifest, write_resolved
from .qc import ensure_qc_provenance, resolve_threshold_rows
from .runtime import assert_config_matches_runtime, print_runtime_identity, record_runtime_provenance
from .slurm import active_cleangene_jobs_for_run, cancel_jobs, sbatch_cmd, submit, user_queue_snapshot
from .task_store import build_isolate_task_store
from .tools import ToolResolutionError, executable_version, resolve_checkm2_executable
from .util import atomic_json, command_exists, load_json, read_tsv, safe_name, sha256, write_tsv
from .utils_cli import add_utils_parser
from .ux import clean_gene_banner, submitted, spinner, waiting
from .workers import cleanup_trimmed_fastqs, compress_completed_outputs, dispatch, global_preflight, invalidate_legacy_identity_metrics as _invalidate_legacy_identity_metrics, invalidate_legacy_isolate_qc as _invalidate_legacy_isolate_qc, needs_checkm2, needs_kraken, run_resume_maintenance

def _checkm2_limited_env() -> dict[str,str]:
    env=os.environ.copy()
    for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","TF_NUM_INTRAOP_THREADS","TF_NUM_INTEROP_THREADS"):
        env[key]="1"
    return env

class LauncherTiming:
    def __init__(self):
        self.run: Path | None = None
        self.rows: list[list[str]] = []
        self.started = time.monotonic()
    def set_run(self, run: Path) -> None:
        self.run = run
    def timed(self, phase: str, func):
        started = time.monotonic()
        try: return func()
        finally: self.rows.append([phase, f"{time.monotonic() - started:.6f}"])
    def write(self) -> None:
        if not self.run: return
        self.rows.append(["total", f"{time.monotonic() - self.started:.6f}"])
        write_tsv(self.run/"logs"/"launcher_timing.tsv",["phase","seconds"],self.rows)

def apply_cli_overrides(cfg: dict[str,str], args) -> dict[str,str]:
    cfg=dict(cfg)
    assembler=getattr(args,"assembler",None)
    if assembler:
        cfg["ASSEMBLER"]=assembler
        cfg["SKIP_SHOVILL"]="true" if assembler=="off" else "false"
        if assembler in {"spades","off"}:
            cfg["SKIP_TRIM"]="true"
            cfg["READ_TRIMMING_MODE"]="off"
    if getattr(args,"skip_trim",False):
        cfg["SKIP_TRIM"]="true"
        cfg["READ_TRIMMING_MODE"]="off"
    if getattr(args,"skip_shovill",False):
        cfg["SKIP_SHOVILL"]="true"
        cfg["SKIP_TRIM"]="true"
        cfg["READ_TRIMMING_MODE"]="off"
    compress=getattr(args,"compress_assembly_outputs",None)
    if compress: cfg["COMPRESS_ASSEMBLY_OUTPUTS"]=compress
    compress_annotation=getattr(args,"compress_annotation_outputs",None)
    if compress_annotation: cfg["COMPRESS_ANNOTATION_OUTPUTS"]=compress_annotation
    if getattr(args,"cleanup_trimmed_fastq",False): cfg["CLEANUP_TRIMMED_FASTQ"]="true"
    return cfg

def validate_fastq_inputs(rows: list[dict[str,str]]) -> None:
    """Lightweight structural validation only; no remote filesystem access."""
    errors=[]
    for row in rows:
        if row.get("raw_bam"): continue
        if any(row.get(column,"").strip() for column in ("assembly","gff","protein_fasta","checkm2_report","pangenome_dir")): continue
        r1=row.get("R1","").strip(); r2=row.get("R2","").strip()
        if not r1 or not r2:
            errors.append(f"{row['isolate_id']}: missing R1/R2")
    if errors:
        raise SystemExit("Manifest input validation failed:\n" + "\n".join(errors[:20]))

def make_run(manifest: Path, analysis_root: Path, cfg: dict[str,str], run_id: str | None) -> Path:
    rid=run_id or datetime.now().strftime("%y%m%d_%H%M%S_cleangene"); run=analysis_root/"runs"/rid
    (run/"provenance").mkdir(parents=True,exist_ok=True); (run/"state").mkdir(exist_ok=True); (run/"logs"/"slurm").mkdir(parents=True,exist_ok=True)
    rows=load_manifest(manifest); validate_fastq_inputs(rows); write_resolved(run/"provenance"/"manifest.tsv",rows); atomic_json(run/"provenance"/"resolved_config.json",cfg); atomic_json(run/"provenance"/"inputs.json",{"manifest":str(manifest.resolve()),"manifest_sha256":sha256(manifest),"created":datetime.now().isoformat()})
    write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],([r["group_id"],r["isolate_id"]] for r in rows))
    counts={}
    for r in rows: counts[r["group_id"]]=counts.get(r["group_id"],0)+1
    write_tsv(run/"state"/"group_tasks.tsv",["group_id","n_isolates","group_size_class"],([g,counts[g],"unresolved"] for g in groups(rows)))
    return run

def load_existing(root: Path, run_id: str) -> Path:
    run=root/"runs"/run_id
    if not (run/"provenance"/"resolved_config.json").is_file(): raise SystemExit(f"Run not found: {run}")
    return run

def validate_run_dir(run: Path) -> None:
    required=[run/"provenance"/"resolved_config.json",run/"provenance"/"manifest.tsv",run/"state"/"isolate_tasks.tsv"]
    missing=[str(p) for p in required if not p.is_file()]
    if missing: raise SystemExit("Cannot resume run; missing metadata:\n" + "\n".join(missing))

def refresh_resume_config(run: Path, config: Path | None) -> dict[str,str]:
    current=load_json(run/"provenance"/"resolved_config.json")
    if not config: return current
    updated=read_env(config)
    if current.get("KRAKEN2_DB") and not updated.get("KRAKEN2_DB"): updated["KRAKEN2_DB"]=current["KRAKEN2_DB"]
    if current.get("CHECKM2_DB") and not updated.get("CHECKM2_DB"): updated["CHECKM2_DB"]=current["CHECKM2_DB"]
    if current.get("CHECKM2_EXECUTABLE") and not updated.get("CHECKM2_EXECUTABLE"): updated["CHECKM2_EXECUTABLE"]=current["CHECKM2_EXECUTABLE"]
    if current.get("CHECKM2_VERSION") and not updated.get("CHECKM2_VERSION"): updated["CHECKM2_VERSION"]=current["CHECKM2_VERSION"]
    if current.get("QC_PROFILE_FILE"): updated["QC_PROFILE_FILE"]=current["QC_PROFILE_FILE"]
    backup=run/"provenance"/"resolved_config.pre_resume.json"
    if not backup.is_file(): shutil.copy2(run/"provenance"/"resolved_config.json",backup)
    atomic_json(run/"provenance"/"resolved_config.json",updated)
    return updated

def invalidate_legacy_identity_metrics(run: Path, cfg: dict[str,str]) -> int:
    return _invalidate_legacy_identity_metrics(run,cfg)

def invalidate_legacy_isolate_qc(run: Path) -> int:
    return _invalidate_legacy_isolate_qc(run)

def latest_run(root: Path) -> Path:
    runs=sorted((root/"runs").glob("*"),key=lambda p:p.stat().st_mtime,reverse=True)
    if not runs: raise SystemExit(f"No runs found under {root/'runs'}")
    return runs[0]

def check(args) -> int:
    rows=load_manifest(args.manifest); cfg=apply_cli_overrides(read_env(args.config),args); required=[]
    profile=Path(cfg["QC_PROFILE_FILE"]).expanduser().resolve() if cfg.get("QC_PROFILE_FILE","").strip() else None
    resolve_threshold_rows(rows,cfg,profile)
    assembler=assembler_mode(cfg)
    if assembler=="off":
        if any(r.get("raw_bam") for r in rows): required.append("samtools")
    else:
        required=["shovill" if assembler=="shovill" else "spades.py","prokka","panaroo","bwa","samtools","bcftools","minimap2"]
    needs_kraken=cfg["TAXONOMY_MODE"] not in {"off","auto"} or any(r.get("grouping_source")=="kraken_pending" for r in rows)
    if needs_kraken: required.append("kraken2")
    mode=checkm2_mode(cfg)
    trim_mode="off" if assembler in {"spades","off"} or truthy(cfg.get("SKIP_TRIM","false")) else cfg.get("READ_TRIMMING_MODE","auto")
    if trim_mode=="always": required.append("fastp")
    missing=[x for x in required if not command_exists(x)]
    if mode=="required":
        try: resolve_checkm2_executable(cfg.get("CHECKM2_EXECUTABLE",""))
        except ToolResolutionError: missing.append("checkm2")
    checkm2_db=Path(cfg.get("CHECKM2_DB","")).expanduser() if cfg.get("CHECKM2_DB","").strip() else None
    if mode=="required" and checkm2_db and not checkm2_db.is_file(): missing.append("CHECKM2_DB")
    print(f"manifest: {len(rows)} isolates / {len(groups(rows))} groups")
    print(f"inputs: {sum(1 for r in rows if r.get('raw_bam'))} paired uBAM / {sum(1 for r in rows if r.get('R1') and r.get('R2'))} FASTQ-pair rows")
    print("tools: " + ("OK" if not missing else "missing " + ", ".join(missing)))
    if trim_mode=="auto" and not command_exists("fastp"): print("warning: fastp not found; adapter trimming check will be skipped")
    if needs_kraken and not cfg.get("KRAKEN2_DB"): print("warning: KRAKEN2_DB is not configured; run will use KRAKEN2_AUTO_DOWNLOAD if enabled")
    return 0 if not missing else 2

def preflight_runtime(run: Path, cfg: dict[str,str]) -> dict[str,str]:
    rows=read_tsv(run/"provenance"/"manifest.tsv")
    if not needs_checkm2(rows,cfg):
        return cfg
    try:
        executable=resolve_checkm2_executable(cfg.get("CHECKM2_EXECUTABLE",""))
        version=executable_version(executable,"CheckM2")
    except ToolResolutionError as error:
        raise SystemExit("CheckM2 is missing from the CleanGene runtime environment before controller submission:\n" + str(error))
    updated={**cfg,"CHECKM2_EXECUTABLE":str(executable),"CHECKM2_VERSION":version}
    atomic_json(run/"provenance"/"resolved_config.json",updated)
    return updated

def _destination_is_writable(path: Path) -> bool:
    candidate=path
    while not candidate.exists() and candidate != candidate.parent:
        candidate=candidate.parent
    return candidate.is_dir() and os.access(candidate,os.W_OK)

def _doctor_config_errors(cfg: dict[str,str]) -> list[str]:
    errors=[]
    try: assembler_mode(cfg)
    except SystemExit as error: errors.append(str(error))
    try: checkm2_mode(cfg)
    except SystemExit as error: errors.append(str(error))
    enums={
        "TAXONOMY_MODE":{"auto","identify","contamination","kraken2","off"},
        "READ_TRIMMING_MODE":{"off","auto","always"},
        "KRAKEN2_DB_ACCESS":{"auto","copy","mmap","direct"},
    }
    for key,allowed in enums.items():
        value=cfg.get(key,"").strip().lower()
        if value not in allowed: errors.append(f"{key} must be one of {', '.join(sorted(allowed))}; got {value!r}")
    positive=("CPUS","SLURM_CPUS","SLURM_ARRAY_CHUNK_SIZE","SLURM_MAX_OUTSTANDING_CHUNKS","SLURM_USER_JOB_LIMIT","SLURM_CONTROLLER_CPUS","KRAKEN2_DB_CPUS","CHECKM2_CPUS","CHECKM2_PREDICT_CPUS","CHECKM2_MAX_INFLIGHT","PREFLIGHT_FILE_CHECK_WORKERS","VALIDATION_CPUS")
    nonnegative=("SLURM_MAX_PARALLEL","SLURM_PREPROCESS_MAX_INFLIGHT","SLURM_VALIDATION_MAX_INFLIGHT","SLURM_GROUP_MAX_INFLIGHT","SLURM_JOB_HEADROOM","SLURM_USER_CPU_LIMIT","SLURM_CPU_HEADROOM","SLURM_POLL_SECONDS")
    for key in positive+nonnegative:
        try: value=int(cfg.get(key,""))
        except ValueError: errors.append(f"{key} must be an integer; got {cfg.get(key,'')!r}"); continue
        if key in positive and value<=0: errors.append(f"{key} must be greater than zero; got {value}")
        if key in nonnegative and value<0: errors.append(f"{key} must be zero or greater; got {value}")
    for key,value in cfg.items():
        if key.endswith("_MEM") and value and not re.fullmatch(r"[1-9][0-9]*(?:\.[0-9]+)?[KMGT]",value,re.IGNORECASE):
            errors.append(f"{key} must be a positive Slurm memory value such as 32G; got {value!r}")
        if key.endswith("_TIME") and value and not re.fullmatch(r"(?:[0-9]+-)?[0-9]+:[0-5][0-9]:[0-5][0-9]",value):
            errors.append(f"{key} must use Slurm D-HH:MM:SS or HH:MM:SS format; got {value!r}")
    try:
        if int(cfg["SLURM_JOB_HEADROOM"])>=int(cfg["SLURM_USER_JOB_LIMIT"]): errors.append("SLURM_JOB_HEADROOM must be smaller than SLURM_USER_JOB_LIMIT")
    except (KeyError,ValueError): pass
    return errors

def doctor(args) -> int:
    cfg=apply_cli_overrides(read_env(args.config),args)
    assert_config_matches_runtime(args.config,cfg)
    print_runtime_identity(args.config)
    failures=0
    prefix_value=os.environ.get("CONDA_PREFIX","")
    prefix=Path(prefix_value).expanduser().resolve() if prefix_value else None
    env_name=os.environ.get("CONDA_DEFAULT_ENV","") or (prefix.name if prefix else "")
    python=Path(sys.executable).resolve()
    import cleangene
    package_path=Path(cleangene.__file__).resolve().parent
    checkout=Path(__file__).resolve().parents[2]
    print(f"CleanGene source: READY {checkout}")
    print(f"CleanGene commit: {__import__('cleangene.runtime',fromlist=['cleangene_commit']).cleangene_commit(checkout)}")
    if package_path == checkout/"src"/"cleangene":
        print("CleanGene editable/runtime match: READY")
    else:
        failures+=1; print(f"CleanGene editable/runtime match: ERROR package={package_path} expected={checkout/'src'/'cleangene'}")
    if env_name=="cleangene" and prefix is not None and prefix in python.parents: print("CleanGene Python environment: OK")
    else:
        failures+=1; print(f"CleanGene Python environment: ERROR expected cleangene Python under its active Conda prefix, found environment={env_name or '<unknown>'} python={python}. Fix: conda activate cleangene")
    config_errors=_doctor_config_errors(cfg)
    if config_errors:
        failures+=len(config_errors)
        for error in config_errors: print(f"Configuration: ERROR {error}")
        target=args.config or Path("config/cleangene.arc.local.env")
        print(f"Configuration fix: edit {target}, then run mamba run -n cleangene cleangene doctor --config {target}")
    else: print("Configuration: OK")
    for tool in ("shovill","spades.py","prokka","panaroo","bwa","samtools","bcftools","minimap2","fastp","kraken2"):
        if command_exists(tool): print(f"Primary tool {tool}: OK")
        else:
            failures+=1; print(f"Primary tool {tool}: ERROR missing. Fix: bash scripts/install_or_update.sh --recreate")
    if cfg.get("CHECKM2_MODE","required").strip().lower()=="required":
        executable=None; version=""
        try:
            executable=resolve_checkm2_executable(cfg.get("CHECKM2_EXECUTABLE",""))
            version=executable_version(executable,"CheckM2")
            print(f"CheckM2 executable: READY {executable}")
            print(f"CheckM2 version: {version}")
            capabilities=checkm2_predict_capabilities(executable)
            print("CheckM2 predict CLI: READY")
            print(f"CheckM2 cleanup option: {capabilities.cleanup_option}")
        except ToolResolutionError as error:
            failures+=1; print(f"CheckM2 executable: ERROR {error}. Fix: bash scripts/install_or_update.sh --recreate")
        except CheckM2DbError as error:
            failures+=1; print(f"CheckM2 predict CLI: ERROR {error}. Fix: bash scripts/install_or_update.sh --recreate")
        try:
            resolution=resolve_checkm2_db(cfg,allow_download=False)
            print(f"CheckM2 database: READY {resolution.path}")
            if getattr(args,"deep_checkm2",False) and executable is not None:
                with tempfile.TemporaryDirectory() as tmp:
                    tmpdir=Path(tmp); env=_checkm2_limited_env()
                    testrun=subprocess.run(checkm2_testrun_command(executable,resolution.path,1),cwd=tmpdir,capture_output=True,text=True,env=env)
                    if testrun.returncode:
                        raise CheckM2DbError((testrun.stderr or testrun.stdout or "").strip() or f"testrun exited {testrun.returncode}")
                    print("CheckM2 testrun: READY")
                    genome=bundled_test_genome(executable)
                    input_path=checkm2_named_input_link(genome,tmpdir/"input","cleangene_checkm2_smoke")
                    out=tmpdir/"predict"
                    command=checkm2_predict_command(executable,input_path,out,resolution.path,1,checkm2_predict_capabilities(executable))
                    predict=subprocess.run(command,cwd=tmpdir,capture_output=True,text=True,env=env)
                    if predict.returncode:
                        raise CheckM2DbError((predict.stderr or predict.stdout or "").strip() or f"predict exited {predict.returncode}")
                    parse_checkm2_quality_report(out/"quality_report.tsv","cleangene_checkm2_smoke")
                    print("CheckM2 production predict smoke test: READY")
        except CheckM2DbNotReady:
            root=checkm2_database_root(cfg)
            if not truthy(cfg.get("CHECKM2_AUTO_DOWNLOAD","true")):
                failures+=1; print(f"CheckM2 database: ERROR not present and CHECKM2_AUTO_DOWNLOAD=false. Fix: set CHECKM2_AUTO_DOWNLOAD=true in {args.config or 'the config'}")
            elif _destination_is_writable(root):
                print(f"CheckM2 database: not present; will be created by checkm2_db_setup | root={root}")
            else:
                failures+=1; print(f"CheckM2 database: ERROR managed root is not writable: {root}. Fix: set CLEANGENE_DATABASE_ROOT in {args.config or 'the config'} to a writable shared directory")
        except CheckM2DbError as error:
            failures+=1; print(f"CheckM2 database: ERROR {error}")
    if cfg.get("TAXONOMY_MODE","auto")!="off":
        try:
            from .kraken import Kraken2DbNotReady, managed_kraken2_db_path, resolve_kraken2_db
            resolution=resolve_kraken2_db(cfg,allow_download=False)
            print(f"Kraken2 configuration: OK database={resolution.path}")
        except Kraken2DbNotReady as error:
            _,_,root=managed_kraken2_db_path(cfg)
            if truthy(cfg.get("KRAKEN2_AUTO_DOWNLOAD","true")) and _destination_is_writable(root): print(f"Kraken2 configuration: OK managed database will be created by kraken_db_setup | root={root}")
            else:
                failures+=1; print(f"Kraken2 configuration: ERROR {error}. Fix: set CLEANGENE_DATABASE_ROOT in {args.config or 'the config'} to a writable shared directory")
        except Exception as error:
            failures+=1; print(f"Kraken2 configuration: ERROR {error}")
    if command_exists("sbatch"): print("SLURM: READY")
    else:
        failures+=1; print("Slurm sbatch: ERROR not found. Fix: run setup and CleanGene from an ARC login node with Slurm commands available")
    return 0 if failures==0 else 2

def estimate(args) -> int:
    rows=load_manifest(args.manifest); total=sum(Path(r[k]).stat().st_size for r in rows for k in (("raw_bam",) if r.get("raw_bam") else ("R1","R2"))); ng=len(groups(rows)); n=len(rows)
    print(json.dumps({"isolates":n,"groups":ng,"compressed_read_bytes":total,"slurm_preprocess_tasks":n,"slurm_panaroo_tasks":ng,"slurm_validation_tasks":n,"slurm_reduce_tasks":ng},indent=2)); return 0

def local(run: Path) -> None:
    cfg,rows={**DEFAULTS,**load_json(run/"provenance"/"resolved_config.json")},read_tsv(run/"provenance"/"manifest.tsv")
    if needs_kraken(rows,cfg): dispatch("kraken_db_setup",run,None)
    cfg,rows={**DEFAULTS,**load_json(run/"provenance"/"resolved_config.json")},read_tsv(run/"provenance"/"manifest.tsv")
    if needs_checkm2(rows,cfg): dispatch("checkm2_db_setup",run,None)
    ni=len((run/"state"/"isolate_tasks.tsv").read_text().splitlines())-1
    for i in range(ni): dispatch("preprocess",run,i)
    dispatch("resolve_groups",run,None)
    ng=len((run/"state"/"group_tasks.tsv").read_text().splitlines())-1
    for i in range(ng): dispatch("panaroo",run,i)
    for i in range(ng): dispatch("prepare_validation",run,i)
    for i in range(ni): dispatch("validate",run,i)
    for i in range(ng): dispatch("reduce",run,i)
    for i in range(ng): dispatch("plot",run,i)
    dispatch("summary",run,None)

def slurm(run: Path, cfg: dict[str,str], dry: bool) -> str:
    exe=f"{shlex_quote(sys.executable)} -m cleangene _worker"
    base=dict(account=cfg["SLURM_ACCOUNT"],partition=cfg["SLURM_PARTITION"])
    wrap=f"{exe} --stage slurm_controller --run-dir {shlex_quote(str(run))} --index 0"
    log=run/"logs"/"slurm"/"controller.%j.log"
    cmd=sbatch_cmd(name="cg-controller",wrap=wrap,cpus=cfg["SLURM_CONTROLLER_CPUS"],mem=cfg["SLURM_CONTROLLER_MEM"],time=cfg["SLURM_CONTROLLER_TIME"],log=log,**base)
    return submit(cmd,dry)

def _guard_resume_active_jobs(run: Path, *, cancel_active: bool = False) -> None:
    active=active_cleangene_jobs_for_run(run)
    if not active:
        return
    job_ids=sorted({str(entry["job_id"]) for entry in active})
    command="scancel " + " ".join(job_ids)
    if cancel_active:
        cancel_jobs(job_ids)
        print(waiting(f"Canceled active CleanGene jobs for this run: {' '.join(job_ids)}"),flush=True)
        return
    raise SystemExit(
        "Refusing to resume while active CleanGene jobs reference this run directory.\n"
        f"Run directory: {run}\n"
        f"Active job IDs: {', '.join(job_ids)}\n"
        f"Cancel only these jobs with:\n  {command}\n"
        "Then rerun the resume command."
    )

def shlex_quote(x:str)->str:
    import shlex; return shlex.quote(x)

def run_command(args) -> int:
    timing=LauncherTiming()
    print(clean_gene_banner())
    print(submitted("Welcome to CleanGene, Your Grace."))
    cfg=timing.timed("parse_config",lambda: apply_cli_overrides(read_env(args.config),args)); assert_config_matches_runtime(args.config,cfg); root=args.analysis_root.expanduser().resolve(); root.mkdir(parents=True,exist_ok=True)
    if args.resume:
        run=load_existing(root,args.resume); timing.set_run(run); cfg=apply_cli_overrides(refresh_resume_config(run,args.config),args)
        if args.skip_trim or args.skip_shovill or getattr(args,"assembler",None) or args.compress_assembly_outputs or args.compress_annotation_outputs or getattr(args,"cleanup_trimmed_fastq",False): atomic_json(run/"provenance"/"resolved_config.json",cfg)
    else:
        run_id=args.run_id or datetime.now().strftime("%y%m%d_%H%M%S_cleangene"); run=root/"runs"/run_id
    print(f"Run directory: {run}")
    with spinner("Getting ready to submit"):
        if not args.resume:
            run=timing.timed("create_run",lambda: make_run(args.manifest,root,cfg,run_id)); timing.set_run(run)
        else: print(waiting("step=resume: submitting controller; legacy checks will run inside the controller job"),flush=True)
        cfg={**DEFAULTS,**load_json(run/"provenance"/"resolved_config.json")}
        record_runtime_provenance(run,cfg)
        if args.profile=="local":
            global_preflight(run)
            if args.resume: run_resume_maintenance(run,cfg)
            local(run)
        else:
            if args.resume: _guard_resume_active_jobs(run,cancel_active=getattr(args,"cancel_active",False))
            controller_job_id=timing.timed("submit_controller",lambda: slurm(run,cfg,args.dry_run))
    timing.write()
    if args.profile=="slurm":
        print(submitted(f"CleanGene run created: {run.name}"))
        print(submitted(f"Controller submitted: {controller_job_id}"))
        print(submitted(f"Run directory: {run}"))
    return 0

def resume_command(args) -> int:
    timing=LauncherTiming()
    print(clean_gene_banner())
    print(submitted("Welcome to CleanGene, Your Grace."))
    if args.run_dir: run=args.run_dir.expanduser().resolve()
    else:
        root=(args.analysis_root or Path.cwd()).expanduser().resolve()
        run=latest_run(root) if args.latest else load_existing(root,args.run)
    validate_run_dir(run)
    timing.set_run(run)
    cfg=apply_cli_overrides(refresh_resume_config(run,args.config),args)
    assert_config_matches_runtime(args.config,cfg)
    if args.skip_trim or args.skip_shovill or getattr(args,"assembler",None) or args.compress_assembly_outputs or args.compress_annotation_outputs or getattr(args,"cleanup_trimmed_fastq",False): atomic_json(run/"provenance"/"resolved_config.json",cfg)
    cfg={**DEFAULTS,**load_json(run/"provenance"/"resolved_config.json")}
    record_runtime_provenance(run,cfg)
    print(f"Run directory: {run}")
    with spinner("Getting ready to submit"):
        print(waiting("step=resume: submitting controller; legacy checks will run inside the controller job"),flush=True)
        _guard_resume_active_jobs(run,cancel_active=getattr(args,"cancel_active",False))
        controller_job_id=timing.timed("submit_controller",lambda: slurm(run,cfg,args.dry_run))
    timing.write()
    print(submitted(f"CleanGene run created: {run.name}"))
    print(submitted(f"Controller submitted: {controller_job_id}"))
    print(submitted(f"Run directory: {run}"))
    return 0

def cleanup_command(args) -> int:
    run=args.run_dir.expanduser().resolve(); validate_run_dir(run)
    if not (run/"state"/"summary.done.json").is_file():
        raise SystemExit("Cleanup refused: the run has not finished (state/summary.done.json is missing)")
    result=cleanup_trimmed_fastqs(run,dry_run=args.dry_run)
    action="would reclaim" if args.dry_run else "reclaimed"
    print(f"{action} {int(result['bytes_reclaimed'])/1024**3:.2f} GiB")
    for status,count in result["counts"].items(): print(f"{status}: {count}")
    if args.dry_run: print("No files were changed. Re-run without --dry-run to apply cleanup.")
    else: print(f"Report: {run/'results'/'cohort'/'fastq_cleanup.tsv'}")
    return 0

def reconcile_preprocess_command(args) -> int:
    run=args.run_dir.expanduser().resolve(); validate_run_dir(run)
    apply_changes=bool(args.apply)
    cfg={**DEFAULTS,**load_json(run/"provenance"/"resolved_config.json")}
    rows=read_tsv(run/"state"/"isolate_tasks.tsv")
    try: snapshot=user_queue_snapshot()
    except (OSError, subprocess.SubprocessError): snapshot={"total":0,"jobs":{},"entries":[]}
    counts=reconcile_preprocess_outputs(run,cfg,rows,snapshot,apply=apply_changes)
    if args.compress_safe:
        if not apply_changes: raise SystemExit("--compress-safe requires --apply")
        compressed=compress_completed_outputs(run)
        print(f"compressed_outputs={compressed.get('files_compressed',0)} bytes_reclaimed={compressed.get('bytes_reclaimed',0)}")
    print("preprocess_reconciliation: " + " ".join(f"{key}={counts[key]}" for key in ("total","marker_complete","output_recovered","active","incomplete","inconsistent")))
    print(f"Report: {run/'state'/'preprocess_reconciliation.tsv'}")
    if not apply_changes: print("Dry run only. Re-run with --apply to recreate safe missing markers or quarantine stale markers.")
    if args.require_all and (counts.get("incomplete",0) or counts.get("inconsistent",0)): return 1
    return 0

def exclude_command(args) -> int:
    run=args.run_dir.expanduser().resolve(); validate_run_dir(run)
    downstream=[]
    for stage in ("panaroo","prepare_validation","validate","reduce","plot"):
        downstream.extend((run/"state"/stage).glob("*.done.json"))
    if downstream or (run/"state"/"summary.done.json").is_file() or any((run/"results"/"groups").glob("*/02_pangenome/panaroo/gene_presence_absence.csv")):
        raise SystemExit("Exclusion refused: downstream pangenome work has already started. Create a new run from a filtered manifest to avoid mixing old and new cohort results.")
    requested=list(args.samples or [])
    if args.samples_file:
        values=[line.strip().split("\t",1)[0] for line in args.samples_file.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
        if values and values[0].lower() in {"isolate_id","sample_id"}: values=values[1:]
        requested.extend(values)
    requested=list(dict.fromkeys(requested))
    if not requested: raise SystemExit("Provide at least one isolate with --samples or --samples-file")
    manifest=run/"provenance"/"manifest.tsv"; rows=read_tsv(manifest); known={row["isolate_id"] for row in rows}; missing=[sample for sample in requested if sample not in known]
    if missing: raise SystemExit("Isolates not found in run manifest: " + ", ".join(missing[:20]))
    selected=set(requested)
    for row in rows:
        if row["isolate_id"] in selected:
            row["user_excluded"]="true"
            marker=run/"state"/"preprocess"/f"{safe_name(row['isolate_id'])}.done.json"
            if marker.is_file(): marker.unlink()
    write_resolved(manifest,rows)
    build_isolate_task_store(run,rows)
    write_tsv(run/"provenance"/"user_exclusions.tsv",["isolate_id","status"],([sample,"user_excluded"] for sample in requested))
    print(f"Marked {len(requested)} isolates as user-excluded without changing task indices.")
    print(f"Report: {run/'provenance'/'user_exclusions.tsv'}")
    return 0

def main(argv=None) -> int:
    p=argparse.ArgumentParser(prog="cleangene"); sub=p.add_subparsers(dest="cmd",required=True)
    c=sub.add_parser("check"); c.add_argument("--manifest",type=Path,required=True); c.add_argument("--config",type=Path); c.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); c.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); c.add_argument("--assembler",choices=("shovill","spades","off")); c.add_argument("--compress-assembly-outputs","--compress_assembly_outputs",dest="compress_assembly_outputs",choices=("off","intermediates","all")); c.add_argument("--compress-annotation-outputs","--compress_annotation_outputs",dest="compress_annotation_outputs",choices=("off","nonessential")); c.add_argument("--cleanup-trimmed-fastq","--cleanup_trimmed_fastq",dest="cleanup_trimmed_fastq",action="store_true"); c.set_defaults(func=check)
    e=sub.add_parser("estimate"); e.add_argument("--manifest",type=Path,required=True); e.set_defaults(func=estimate)
    d=sub.add_parser("doctor"); d.add_argument("--config",type=Path); d.add_argument("--manifest",type=Path); d.add_argument("--deep-checkm2",action="store_true"); d.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); d.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); d.add_argument("--assembler",choices=("shovill","spades","off")); d.set_defaults(func=doctor)
    r=sub.add_parser("run"); r.add_argument("--manifest",type=Path); r.add_argument("--analysis-root",type=Path,required=True); r.add_argument("--config",type=Path); r.add_argument("--profile",choices=("local","slurm"),default="slurm"); r.add_argument("--dry-run",action="store_true"); r.add_argument("--run-id"); r.add_argument("--resume"); r.add_argument("--cancel-active",action="store_true"); r.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); r.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); r.add_argument("--assembler",choices=("shovill","spades","off")); r.add_argument("--compress-assembly-outputs","--compress_assembly_outputs",dest="compress_assembly_outputs",choices=("off","intermediates","all")); r.add_argument("--compress-annotation-outputs","--compress_annotation_outputs",dest="compress_annotation_outputs",choices=("off","nonessential")); r.add_argument("--cleanup-trimmed-fastq","--cleanup_trimmed_fastq",dest="cleanup_trimmed_fastq",action="store_true"); r.set_defaults(func=run_command)
    rs=sub.add_parser("resume"); rs.add_argument("--run"); rs.add_argument("--run-dir",type=Path); rs.add_argument("--latest",action="store_true"); rs.add_argument("--analysis-root",type=Path); rs.add_argument("--config",type=Path); rs.add_argument("--dry-run",action="store_true"); rs.add_argument("--cancel-active",action="store_true"); rs.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); rs.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); rs.add_argument("--assembler",choices=("shovill","spades","off")); rs.add_argument("--compress-assembly-outputs","--compress_assembly_outputs",dest="compress_assembly_outputs",choices=("off","intermediates","all")); rs.add_argument("--compress-annotation-outputs","--compress_annotation_outputs",dest="compress_annotation_outputs",choices=("off","nonessential")); rs.add_argument("--cleanup-trimmed-fastq","--cleanup_trimmed_fastq",dest="cleanup_trimmed_fastq",action="store_true"); rs.set_defaults(func=resume_command)
    cl=sub.add_parser("cleanup",help="replace retained trimmed FASTQs with links to original FASTQ inputs"); cl.add_argument("--run-dir",type=Path,required=True); cl.add_argument("--dry-run",action="store_true"); cl.set_defaults(func=cleanup_command)
    rp=sub.add_parser("reconcile-preprocess",help="audit or repair missing preprocess markers from existing qc.tsv outputs"); rp.add_argument("--run-dir",type=Path,required=True); rp.add_argument("--dry-run",action="store_true",default=True); rp.add_argument("--apply",action="store_true"); rp.add_argument("--compress-safe",action="store_true"); rp.add_argument("--require-all",action="store_true"); rp.set_defaults(func=reconcile_preprocess_command)
    x=sub.add_parser("exclude",help="exclude isolates safely before downstream pangenome stages start"); x.add_argument("--run-dir",type=Path,required=True); x.add_argument("--samples",nargs="*"); x.add_argument("--samples-file",type=Path); x.set_defaults(func=exclude_command)
    w=sub.add_parser("_worker"); w.add_argument("--stage",required=True); w.add_argument("--run-dir",type=Path,required=True); w.add_argument("--index",type=int,default=0); w.set_defaults(func=lambda a:(dispatch(a.stage,a.run_dir,a.index),0)[1])
    uw=sub.add_parser("_utils_worker"); uw.add_argument("--request",type=Path,required=True); uw.set_defaults(func=lambda a:(__import__("cleangene.downstream",fromlist=["run_request"]).run_request(a.request),0)[1])
    add_utils_parser(sub)
    args=p.parse_args(argv)
    if args.cmd=="run" and not args.resume and not args.manifest: p.error("run requires --manifest unless --resume is used")
    if args.cmd=="resume" and sum(bool(x) for x in (args.run,args.run_dir,args.latest))!=1: p.error("resume requires exactly one of --run, --run-dir, or --latest")
    if args.cmd=="resume" and args.latest and not args.analysis_root: p.error("resume --latest requires --analysis-root")
    return args.func(args)
