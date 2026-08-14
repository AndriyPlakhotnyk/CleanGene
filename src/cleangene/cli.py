from __future__ import annotations
import argparse, json, os, subprocess, sys
from datetime import datetime
from pathlib import Path
from .config import read_env
from .defaults import SCIENTIFIC_DEFAULTS
from .manifest import groups, load_manifest, write_resolved
from .slurm import sbatch_cmd, submit
from .util import atomic_json, command_exists, safe_name, sha256, write_tsv
from .workers import dispatch

def make_run(manifest: Path, analysis_root: Path, cfg: dict[str,str], run_id: str | None) -> Path:
    rid=run_id or datetime.now().strftime("%y%m%d_%H%M%S_cleangene"); run=analysis_root/"runs"/rid
    (run/"provenance").mkdir(parents=True,exist_ok=True); (run/"state").mkdir(exist_ok=True); (run/"logs"/"slurm").mkdir(parents=True,exist_ok=True)
    rows=load_manifest(manifest); write_resolved(run/"provenance"/"manifest.tsv",rows); atomic_json(run/"provenance"/"resolved_config.json",cfg); atomic_json(run/"provenance"/"inputs.json",{"manifest":str(manifest.resolve()),"manifest_sha256":sha256(manifest),"created":datetime.now().isoformat()})
    write_tsv(run/"state"/"isolate_tasks.tsv",["group_id","isolate_id"],([r["group_id"],r["isolate_id"]] for r in rows)); write_tsv(run/"state"/"group_tasks.tsv",["group_id"],([g] for g in groups(rows)))
    return run

def load_existing(root: Path, run_id: str) -> Path:
    run=root/"runs"/run_id
    if not (run/"provenance"/"resolved_config.json").is_file(): raise SystemExit(f"Run not found: {run}")
    return run

def check(args) -> int:
    rows=load_manifest(args.manifest); cfg=read_env(args.config); required=["shovill","prokka","panaroo","bwa","samtools","bcftools","minimap2"]
    if cfg["TAXONOMY_MODE"]!="off": required.append("kraken2")
    missing=[x for x in required if not command_exists(x)]
    print(f"manifest: {len(rows)} isolates / {len(groups(rows))} groups")
    print("tools: " + ("OK" if not missing else "missing " + ", ".join(missing)))
    if cfg["TAXONOMY_MODE"]!="off" and not cfg.get("KRAKEN2_DB"): print("warning: KRAKEN2_DB must be configured for a real run")
    return 0 if not missing else 2

def estimate(args) -> int:
    rows=load_manifest(args.manifest); total=sum(Path(r[k]).stat().st_size for r in rows for k in ("R1","R2")); ng=len(groups(rows)); n=len(rows)
    print(json.dumps({"isolates":n,"groups":ng,"compressed_read_bytes":total,"slurm_preprocess_tasks":n,"slurm_panaroo_tasks":ng,"slurm_validation_tasks":n,"slurm_reduce_tasks":ng},indent=2)); return 0

def local(run: Path) -> None:
    ni=len((run/"state"/"isolate_tasks.tsv").read_text().splitlines())-1; ng=len((run/"state"/"group_tasks.tsv").read_text().splitlines())-1
    for i in range(ni): dispatch("preprocess",run,i)
    for i in range(ng): dispatch("panaroo",run,i)
    for i in range(ng): dispatch("prepare_validation",run,i)
    for i in range(ni): dispatch("validate",run,i)
    for i in range(ng): dispatch("reduce",run,i)
    dispatch("summary",run,None)

def slurm(run: Path, cfg: dict[str,str], dry: bool) -> None:
    ni=len((run/"state"/"isolate_tasks.tsv").read_text().splitlines())-1; ng=len((run/"state"/"group_tasks.tsv").read_text().splitlines())-1; maxp=cfg["SLURM_MAX_PARALLEL"]; exe=f"{shlex_quote(sys.executable)} -m cleangene _worker"
    base=dict(account=cfg["SLURM_ACCOUNT"],partition=cfg["SLURM_PARTITION"])
    def cmd(stage,array,cpus,mem,time,dep=None):
        idx='${SLURM_ARRAY_TASK_ID}' if array else '0'; wrap=f"{exe} --stage {stage} --run-dir {shlex_quote(str(run))} --index {idx}"; log=run/"logs"/"slurm"/f"{stage}.%A_%a.log"; return sbatch_cmd(name=f"cg-{stage}",wrap=wrap,cpus=cpus,mem=mem,time=time,array=array,dependency=dep,log=log,**base)
    j1=submit(cmd("preprocess",f"0-{ni-1}%{maxp}",cfg["SLURM_CPUS"],cfg["SLURM_MEM"],cfg["SLURM_TIME"]),dry)
    j2=submit(cmd("panaroo",f"0-{ng-1}%{maxp}",cfg["PANAROO_CPUS"],cfg["PANAROO_MEM"],cfg["PANAROO_TIME"],j1),dry)
    j3=submit(cmd("prepare_validation",f"0-{ng-1}%{maxp}",cfg["PANAROO_CPUS"],cfg["PANAROO_MEM"],cfg["PANAROO_TIME"],j2),dry)
    j4=submit(cmd("validate",f"0-{ni-1}%{maxp}",cfg["VALIDATION_CPUS"],cfg["VALIDATION_MEM"],cfg["VALIDATION_TIME"],j3),dry)
    j5=submit(cmd("reduce",f"0-{ng-1}%{maxp}",cfg["SUMMARY_CPUS"],cfg["SUMMARY_MEM"],cfg["SUMMARY_TIME"],j4),dry)
    wrap=f"{exe} --stage summary --run-dir {shlex_quote(str(run))} --index 0"; submit(sbatch_cmd(name="cg-summary",wrap=wrap,cpus=cfg["SUMMARY_CPUS"],mem=cfg["SUMMARY_MEM"],time=cfg["SUMMARY_TIME"],dependency=j5,log=run/"logs"/"slurm"/"summary.%j.log",**base),dry)

def shlex_quote(x:str)->str:
    import shlex; return shlex.quote(x)

def run_command(args) -> int:
    cfg=read_env(args.config); root=args.analysis_root.expanduser().resolve(); root.mkdir(parents=True,exist_ok=True)
    if args.resume: run=load_existing(root,args.resume); cfg=json.load((run/"provenance"/"resolved_config.json").open())
    else: run=make_run(args.manifest,root,cfg,args.run_id)
    print(f"run_dir={run}")
    if args.profile=="local": local(run)
    else: slurm(run,cfg,args.dry_run)
    return 0

def main(argv=None) -> int:
    p=argparse.ArgumentParser(prog="cleangene"); sub=p.add_subparsers(dest="cmd",required=True)
    c=sub.add_parser("check"); c.add_argument("--manifest",type=Path,required=True); c.add_argument("--config",type=Path); c.set_defaults(func=check)
    e=sub.add_parser("estimate"); e.add_argument("--manifest",type=Path,required=True); e.set_defaults(func=estimate)
    r=sub.add_parser("run"); r.add_argument("--manifest",type=Path); r.add_argument("--analysis-root",type=Path,required=True); r.add_argument("--config",type=Path); r.add_argument("--profile",choices=("local","slurm"),default="slurm"); r.add_argument("--dry-run",action="store_true"); r.add_argument("--run-id"); r.add_argument("--resume"); r.set_defaults(func=run_command)
    w=sub.add_parser("_worker"); w.add_argument("--stage",required=True); w.add_argument("--run-dir",type=Path,required=True); w.add_argument("--index",type=int,default=0); w.set_defaults(func=lambda a:(dispatch(a.stage,a.run_dir,a.index),0)[1])
    args=p.parse_args(argv)
    if args.cmd=="run" and not args.resume and not args.manifest: p.error("run requires --manifest unless --resume is used")
    return args.func(args)
