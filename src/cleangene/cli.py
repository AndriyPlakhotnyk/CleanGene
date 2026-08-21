from __future__ import annotations
import argparse, csv, json, os, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path
from .config import read_env
from .defaults import SCIENTIFIC_DEFAULTS
from .manifest import groups, load_manifest, write_resolved
from .slurm import sbatch_cmd, submit
from .util import atomic_json, command_exists, load_json, safe_name, sha256, write_tsv
from .utils_cli import add_utils_parser
from .ux import spinner
from .workers import dispatch

def validate_fastq_inputs(rows: list[dict[str,str]]) -> None:
    """Lightweight structural validation only; no remote filesystem access."""
    errors=[]
    for row in rows:
        if row.get("raw_bam"): continue
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
    backup=run/"provenance"/"resolved_config.pre_resume.json"
    if not backup.is_file(): shutil.copy2(run/"provenance"/"resolved_config.json",backup)
    atomic_json(run/"provenance"/"resolved_config.json",updated)
    return updated

def _float_or_none(value: str | None) -> float | None:
    try:
        return float(value) if value not in {"", "NA", None} else None
    except ValueError:
        return None

def _remove_done(path: Path) -> None:
    if path.is_file(): path.unlink()

def invalidate_legacy_identity_metrics(run: Path, cfg: dict[str,str]) -> int:
    min_breadth=float(cfg["READ_VALIDATION_MIN_BREADTH"]); min_depth=float(cfg["READ_VALIDATION_MIN_MEAN_DEPTH"])
    invalidated=0
    for metrics in (run/"results"/"groups").glob("*/03_read_validation/evidence/*/metrics.tsv"):
        rows=[]
        try:
            with metrics.open(newline="") as handle:
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
        _remove_done(run/"state"/"validate"/f"{iso_safe}.done.json")
        _remove_done(run/"state"/"reduce"/f"{group_safe}.done.json")
        _remove_done(run/"state"/"plot"/f"{group_safe}.done.json")
        _remove_done(run/"state"/"summary.done.json")
        invalidated += 1
    return invalidated

def latest_run(root: Path) -> Path:
    runs=sorted((root/"runs").glob("*"),key=lambda p:p.stat().st_mtime,reverse=True)
    if not runs: raise SystemExit(f"No runs found under {root/'runs'}")
    return runs[0]

def check(args) -> int:
    rows=load_manifest(args.manifest); cfg=read_env(args.config); required=["shovill","prokka","panaroo","bwa","samtools","bcftools","minimap2"]
    needs_kraken=cfg["TAXONOMY_MODE"] not in {"off","auto"} or any(r.get("grouping_source")=="kraken_pending" for r in rows)
    if needs_kraken: required.append("kraken2")
    if cfg.get("READ_TRIMMING_MODE","auto")=="always": required.append("fastp")
    missing=[x for x in required if not command_exists(x)]
    print(f"manifest: {len(rows)} isolates / {len(groups(rows))} groups")
    print(f"inputs: {sum(1 for r in rows if r.get('raw_bam'))} raw BAM / {sum(1 for r in rows if r.get('R1') and r.get('R2'))} FASTQ-pair rows")
    print("tools: " + ("OK" if not missing else "missing " + ", ".join(missing)))
    if cfg.get("READ_TRIMMING_MODE","auto")=="auto" and not command_exists("fastp"): print("warning: fastp not found; adapter trimming check will be skipped")
    if needs_kraken and not cfg.get("KRAKEN2_DB"): print("warning: KRAKEN2_DB is not configured; run will use KRAKEN2_AUTO_DOWNLOAD if enabled")
    return 0 if not missing else 2

def estimate(args) -> int:
    rows=load_manifest(args.manifest); total=sum(Path(r[k]).stat().st_size for r in rows for k in (("raw_bam",) if r.get("raw_bam") else ("R1","R2"))); ng=len(groups(rows)); n=len(rows)
    print(json.dumps({"isolates":n,"groups":ng,"compressed_read_bytes":total,"slurm_preprocess_tasks":n,"slurm_panaroo_tasks":ng,"slurm_validation_tasks":n,"slurm_reduce_tasks":ng},indent=2)); return 0

def local(run: Path) -> None:
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

def shlex_quote(x:str)->str:
    import shlex; return shlex.quote(x)

def run_command(args) -> int:
    print("Welcome to CleanGene")
    cfg=read_env(args.config); root=args.analysis_root.expanduser().resolve(); root.mkdir(parents=True,exist_ok=True)
    if args.resume:
        run=load_existing(root,args.resume); cfg=refresh_resume_config(run,args.config)
    else:
        run_id=args.run_id or datetime.now().strftime("%y%m%d_%H%M%S_cleangene"); run=root/"runs"/run_id
    print(f"Run directory: {run}")
    with spinner("Getting ready to submit"):
        if not args.resume: run=make_run(args.manifest,root,cfg,run_id)
        if args.profile=="local": local(run)
        else: slurm(run,cfg,args.dry_run)
    if args.profile=="slurm": print(f"Run submitted. Please find logs in {run/'logs'/'slurm'}")
    return 0

def resume_command(args) -> int:
    print("Welcome to CleanGene")
    if args.run_dir: run=args.run_dir.expanduser().resolve()
    else:
        root=(args.analysis_root or Path.cwd()).expanduser().resolve()
        run=latest_run(root) if args.latest else load_existing(root,args.run)
    validate_run_dir(run)
    cfg=refresh_resume_config(run,args.config)
    print(f"Run directory: {run}")
    with spinner("Getting ready to submit"):
        invalidated=invalidate_legacy_identity_metrics(run,cfg)
        if invalidated: print(f"invalidated_legacy_identity_metrics={invalidated}")
        slurm(run,cfg,args.dry_run)
    print(f"Run submitted. Please find logs in {run/'logs'/'slurm'}")
    return 0

def main(argv=None) -> int:
    p=argparse.ArgumentParser(prog="cleangene"); sub=p.add_subparsers(dest="cmd",required=True)
    c=sub.add_parser("check"); c.add_argument("--manifest",type=Path,required=True); c.add_argument("--config",type=Path); c.set_defaults(func=check)
    e=sub.add_parser("estimate"); e.add_argument("--manifest",type=Path,required=True); e.set_defaults(func=estimate)
    r=sub.add_parser("run"); r.add_argument("--manifest",type=Path); r.add_argument("--analysis-root",type=Path,required=True); r.add_argument("--config",type=Path); r.add_argument("--profile",choices=("local","slurm"),default="slurm"); r.add_argument("--dry-run",action="store_true"); r.add_argument("--run-id"); r.add_argument("--resume"); r.set_defaults(func=run_command)
    rs=sub.add_parser("resume"); rs.add_argument("--run"); rs.add_argument("--run-dir",type=Path); rs.add_argument("--latest",action="store_true"); rs.add_argument("--analysis-root",type=Path); rs.add_argument("--config",type=Path); rs.add_argument("--dry-run",action="store_true"); rs.set_defaults(func=resume_command)
    w=sub.add_parser("_worker"); w.add_argument("--stage",required=True); w.add_argument("--run-dir",type=Path,required=True); w.add_argument("--index",type=int,default=0); w.set_defaults(func=lambda a:(dispatch(a.stage,a.run_dir,a.index),0)[1])
    uw=sub.add_parser("_utils_worker"); uw.add_argument("--request",type=Path,required=True); uw.set_defaults(func=lambda a:(__import__("cleangene.downstream",fromlist=["run_request"]).run_request(a.request),0)[1])
    add_utils_parser(sub)
    args=p.parse_args(argv)
    if args.cmd=="run" and not args.resume and not args.manifest: p.error("run requires --manifest unless --resume is used")
    if args.cmd=="resume" and sum(bool(x) for x in (args.run,args.run_dir,args.latest))!=1: p.error("resume requires exactly one of --run, --run-dir, or --latest")
    if args.cmd=="resume" and args.latest and not args.analysis_root: p.error("resume --latest requires --analysis-root")
    return args.func(args)
