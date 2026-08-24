from __future__ import annotations
import argparse, json, os, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path
from .config import assembler_mode, checkm2_mode, read_env, truthy
from .defaults import SCIENTIFIC_DEFAULTS
from .manifest import groups, load_manifest, write_resolved
from .qc import ensure_qc_provenance, prepare_qc_provenance, resolve_threshold_rows
from .slurm import sbatch_cmd, submit
from .util import atomic_json, command_exists, load_json, read_tsv, safe_name, sha256, write_tsv
from .utils_cli import add_utils_parser
from .ux import clean_gene_banner, submitted, spinner, waiting, welcome
from .workers import cleanup_trimmed_fastqs, dispatch, invalidate_legacy_identity_metrics as _invalidate_legacy_identity_metrics, invalidate_legacy_isolate_qc as _invalidate_legacy_isolate_qc, run_resume_maintenance

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
        r1=row.get("R1","").strip(); r2=row.get("R2","").strip()
        if not r1 or not r2:
            errors.append(f"{row['isolate_id']}: missing R1/R2")
    if errors:
        raise SystemExit("Manifest input validation failed:\n" + "\n".join(errors[:20]))

def make_run(manifest: Path, analysis_root: Path, cfg: dict[str,str], run_id: str | None) -> Path:
    rid=run_id or datetime.now().strftime("%y%m%d_%H%M%S_cleangene"); run=analysis_root/"runs"/rid
    (run/"provenance").mkdir(parents=True,exist_ok=True); (run/"state").mkdir(exist_ok=True); (run/"logs"/"slurm").mkdir(parents=True,exist_ok=True)
    rows=load_manifest(manifest); validate_fastq_inputs(rows); cfg=prepare_qc_provenance(run,rows,cfg); write_resolved(run/"provenance"/"manifest.tsv",rows); atomic_json(run/"provenance"/"resolved_config.json",cfg); atomic_json(run/"provenance"/"inputs.json",{"manifest":str(manifest.resolve()),"manifest_sha256":sha256(manifest),"created":datetime.now().isoformat()})
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
    if mode=="required": required.append("checkm2")
    trim_mode="off" if assembler in {"spades","off"} or truthy(cfg.get("SKIP_TRIM","false")) else cfg.get("READ_TRIMMING_MODE","auto")
    if trim_mode=="always": required.append("fastp")
    missing=[x for x in required if not command_exists(x)]
    checkm2_db=Path(cfg.get("CHECKM2_DB","")).expanduser() if cfg.get("CHECKM2_DB","").strip() else None
    if mode=="required" and (not checkm2_db or not checkm2_db.is_file()): missing.append("CHECKM2_DB")
    print(f"manifest: {len(rows)} isolates / {len(groups(rows))} groups")
    print(f"inputs: {sum(1 for r in rows if r.get('raw_bam'))} raw BAM / {sum(1 for r in rows if r.get('R1') and r.get('R2'))} FASTQ-pair rows")
    print("tools: " + ("OK" if not missing else "missing " + ", ".join(missing)))
    if trim_mode=="auto" and not command_exists("fastp"): print("warning: fastp not found; adapter trimming check will be skipped")
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
    print(clean_gene_banner())
    print(submitted("Welcome to CleanGene, You Grace."))
    cfg=apply_cli_overrides(read_env(args.config),args); root=args.analysis_root.expanduser().resolve(); root.mkdir(parents=True,exist_ok=True)
    if args.resume:
        run=load_existing(root,args.resume); cfg=apply_cli_overrides(refresh_resume_config(run,args.config),args)
        cfg=ensure_qc_provenance(run,read_tsv(run/"provenance"/"manifest.tsv"),cfg); atomic_json(run/"provenance"/"resolved_config.json",cfg)
        if args.skip_trim or args.skip_shovill or getattr(args,"assembler",None) or args.compress_assembly_outputs or args.compress_annotation_outputs or getattr(args,"cleanup_trimmed_fastq",False): atomic_json(run/"provenance"/"resolved_config.json",cfg)
    else:
        run_id=args.run_id or datetime.now().strftime("%y%m%d_%H%M%S_cleangene"); run=root/"runs"/run_id
    print(f"Run directory: {run}")
    with spinner("Getting ready to submit"):
        if not args.resume: run=make_run(args.manifest,root,cfg,run_id)
        else: print(waiting("step=resume: submitting controller; legacy checks will run inside the controller job"),flush=True)
        if args.profile=="local":
            if args.resume: run_resume_maintenance(run,cfg)
            local(run)
        else: slurm(run,cfg,args.dry_run)
    if args.profile=="slurm":
        print(submitted(f"Run submitted. Please find logs in {run/'logs'/'slurm'}"))
    return 0

def resume_command(args) -> int:
    print(clean_gene_banner())
    print(submitted("Welcome to CleanGene, You Grace."))
    if args.run_dir: run=args.run_dir.expanduser().resolve()
    else:
        root=(args.analysis_root or Path.cwd()).expanduser().resolve()
        run=latest_run(root) if args.latest else load_existing(root,args.run)
    validate_run_dir(run)
    cfg=apply_cli_overrides(refresh_resume_config(run,args.config),args)
    cfg=ensure_qc_provenance(run,read_tsv(run/"provenance"/"manifest.tsv"),cfg); atomic_json(run/"provenance"/"resolved_config.json",cfg)
    if args.skip_trim or args.skip_shovill or getattr(args,"assembler",None) or args.compress_assembly_outputs or args.compress_annotation_outputs or getattr(args,"cleanup_trimmed_fastq",False): atomic_json(run/"provenance"/"resolved_config.json",cfg)
    print(f"Run directory: {run}")
    with spinner("Getting ready to submit"):
        print(waiting("step=resume: submitting controller; legacy checks will run inside the controller job"),flush=True)
        slurm(run,cfg,args.dry_run)
    print(submitted(f"Run submitted. Please find logs in {run/'logs'/'slurm'}"))
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
    write_tsv(run/"provenance"/"user_exclusions.tsv",["isolate_id","status"],([sample,"user_excluded"] for sample in requested))
    print(f"Marked {len(requested)} isolates as user-excluded without changing task indices.")
    print(f"Report: {run/'provenance'/'user_exclusions.tsv'}")
    return 0

def main(argv=None) -> int:
    p=argparse.ArgumentParser(prog="cleangene"); sub=p.add_subparsers(dest="cmd",required=True)
    c=sub.add_parser("check"); c.add_argument("--manifest",type=Path,required=True); c.add_argument("--config",type=Path); c.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); c.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); c.add_argument("--assembler",choices=("shovill","spades","off")); c.add_argument("--compress-assembly-outputs","--compress_assembly_outputs",dest="compress_assembly_outputs",choices=("off","intermediates","all")); c.add_argument("--compress-annotation-outputs","--compress_annotation_outputs",dest="compress_annotation_outputs",choices=("off","nonessential")); c.add_argument("--cleanup-trimmed-fastq","--cleanup_trimmed_fastq",dest="cleanup_trimmed_fastq",action="store_true"); c.set_defaults(func=check)
    e=sub.add_parser("estimate"); e.add_argument("--manifest",type=Path,required=True); e.set_defaults(func=estimate)
    r=sub.add_parser("run"); r.add_argument("--manifest",type=Path); r.add_argument("--analysis-root",type=Path,required=True); r.add_argument("--config",type=Path); r.add_argument("--profile",choices=("local","slurm"),default="slurm"); r.add_argument("--dry-run",action="store_true"); r.add_argument("--run-id"); r.add_argument("--resume"); r.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); r.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); r.add_argument("--assembler",choices=("shovill","spades","off")); r.add_argument("--compress-assembly-outputs","--compress_assembly_outputs",dest="compress_assembly_outputs",choices=("off","intermediates","all")); r.add_argument("--compress-annotation-outputs","--compress_annotation_outputs",dest="compress_annotation_outputs",choices=("off","nonessential")); r.add_argument("--cleanup-trimmed-fastq","--cleanup_trimmed_fastq",dest="cleanup_trimmed_fastq",action="store_true"); r.set_defaults(func=run_command)
    rs=sub.add_parser("resume"); rs.add_argument("--run"); rs.add_argument("--run-dir",type=Path); rs.add_argument("--latest",action="store_true"); rs.add_argument("--analysis-root",type=Path); rs.add_argument("--config",type=Path); rs.add_argument("--dry-run",action="store_true"); rs.add_argument("--skip-trim","--skip_trim",dest="skip_trim",action="store_true"); rs.add_argument("--skip-shovill","--skip_shovill",dest="skip_shovill",action="store_true"); rs.add_argument("--assembler",choices=("shovill","spades","off")); rs.add_argument("--compress-assembly-outputs","--compress_assembly_outputs",dest="compress_assembly_outputs",choices=("off","intermediates","all")); rs.add_argument("--compress-annotation-outputs","--compress_annotation_outputs",dest="compress_annotation_outputs",choices=("off","nonessential")); rs.add_argument("--cleanup-trimmed-fastq","--cleanup_trimmed_fastq",dest="cleanup_trimmed_fastq",action="store_true"); rs.set_defaults(func=resume_command)
    cl=sub.add_parser("cleanup",help="replace retained trimmed FASTQs with links to original FASTQ inputs"); cl.add_argument("--run-dir",type=Path,required=True); cl.add_argument("--dry-run",action="store_true"); cl.set_defaults(func=cleanup_command)
    x=sub.add_parser("exclude",help="exclude isolates safely before downstream pangenome stages start"); x.add_argument("--run-dir",type=Path,required=True); x.add_argument("--samples",nargs="*"); x.add_argument("--samples-file",type=Path); x.set_defaults(func=exclude_command)
    w=sub.add_parser("_worker"); w.add_argument("--stage",required=True); w.add_argument("--run-dir",type=Path,required=True); w.add_argument("--index",type=int,default=0); w.set_defaults(func=lambda a:(dispatch(a.stage,a.run_dir,a.index),0)[1])
    uw=sub.add_parser("_utils_worker"); uw.add_argument("--request",type=Path,required=True); uw.set_defaults(func=lambda a:(__import__("cleangene.downstream",fromlist=["run_request"]).run_request(a.request),0)[1])
    add_utils_parser(sub)
    args=p.parse_args(argv)
    if args.cmd=="run" and not args.resume and not args.manifest: p.error("run requires --manifest unless --resume is used")
    if args.cmd=="resume" and sum(bool(x) for x in (args.run,args.run_dir,args.latest))!=1: p.error("resume requires exactly one of --run, --run-dir, or --latest")
    if args.cmd=="resume" and args.latest and not args.analysis_root: p.error("resume --latest requires --analysis-root")
    return args.func(args)
