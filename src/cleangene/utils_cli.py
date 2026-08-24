from __future__ import annotations
import argparse, json, sys
from datetime import datetime
from pathlib import Path
from .defaults import DEFAULTS
from .downstream import read_ids, read_manifest_rows, resolve_organism, validate_matrix_selection
from .slurm import sbatch_cmd, submit
from .util import atomic_json, read_tsv, safe_name
from .ux import submitted, spinner, welcome

def _run_args(parser: argparse.ArgumentParser) -> None:
    group=parser.add_mutually_exclusive_group(required=True); group.add_argument("--run"); group.add_argument("--run-dir",type=Path); group.add_argument("--latest",action="store_true")
    parser.add_argument("--analysis-root",type=Path); parser.add_argument("--analysis-name"); parser.add_argument("--dry-run",action="store_true")

def _organism_args(parser: argparse.ArgumentParser, required: bool = False) -> None:
    parser.add_argument("--organism",required=required); parser.add_argument("--samples",nargs="*",default=[]); parser.add_argument("--sample-file",type=Path)

def add_utils_parser(sub) -> None:
    root=sub.add_parser("utils",help="submit downstream analyses to SLURM"); utilities=root.add_subparsers(dest="utility",required=True)
    samples=utilities.add_parser("get-samples",aliases=["get_samples"]); _run_args(samples); _organism_args(samples); samples.add_argument("--genes",nargs="+",required=True); samples.add_argument("--status",choices=("present","absent","both"),default="both"); samples.add_argument("--match",choices=("all","any"),default="all")
    diff=utilities.add_parser("get-differential-genes",aliases=["get_differential_genes"]); _run_args(diff); _organism_args(diff); diff.add_argument("--cohort-a",nargs="*",default=[]); diff.add_argument("--cohort-a-file",type=Path); diff.add_argument("--cohort-b",nargs="*",default=[]); diff.add_argument("--cohort-b-file",type=Path); diff.add_argument("--manifest",type=Path); diff.add_argument("--group-column",default="group_id"); diff.add_argument("--group-a-label"); diff.add_argument("--group-b-label"); diff.add_argument("--max-q-value",type=float,default=0.05); diff.add_argument("--min-prevalence-difference",type=float,default=0.10); diff.add_argument("--top",type=int,default=50)
    operon=utilities.add_parser("get-operon",aliases=["get_operon"]); _run_args(operon); _organism_args(operon); operon.add_argument("--genes",nargs="+"); operon.add_argument("--operon-name"); operon.add_argument("--manifest",type=Path)
    variants=utilities.add_parser("get-variants",aliases=["get_variants"]); _run_args(variants); _organism_args(variants); variants.add_argument("--genes",nargs="+"); variants.add_argument("--operon",type=Path); variants.add_argument("--min-similarity",type=float,default=95.0); variants.add_argument("--flanking-genes",type=int,default=0); variants.add_argument("--analyze-flanks",action="store_true")
    diagnostic=utilities.add_parser("diagnose-call",aliases=["diagnose_call"]); _run_args(diagnostic); _organism_args(diagnostic); diagnostic.add_argument("--genes",nargs="+",required=True); diagnostic.add_argument("--max-samples",type=int,default=20); diagnostic.add_argument("--skip-assembly-replay",action="store_true")
    itol=utilities.add_parser("itol",aliases=["get_itol"]); _run_args(itol); _organism_args(itol); itol.add_argument("--genes",nargs="*",default=[]); itol.add_argument("--operon",type=Path); itol.add_argument("--variants",type=Path); itol.add_argument("--color-scheme",choices=("classic","muted","custom"),default="classic"); itol.add_argument("--custom-colors",type=Path)
    root.set_defaults(func=utils_command)

def locate_run(args) -> Path:
    if args.run_dir: run=args.run_dir.expanduser().resolve()
    else:
        if not args.analysis_root: raise SystemExit("--analysis-root is required with --run or --latest")
        runs=args.analysis_root.expanduser().resolve()/"runs"
        if args.latest:
            found=sorted(runs.glob("*"),key=lambda p:p.stat().st_mtime,reverse=True)
            if not found: raise SystemExit(f"No runs found under {runs}")
            run=found[0]
        else: run=runs/args.run
    required=[run/"provenance"/"resolved_config.json",run/"provenance"/"manifest.tsv",run/"state"/"group_tasks.tsv"]
    missing=[str(x) for x in required if not x.is_file()]
    if missing: raise SystemExit("CleanGene run metadata is missing:\n"+"\n".join(missing))
    return run

def _samples(args) -> list[str]:
    return list(dict.fromkeys([*getattr(args,"samples",[]),*read_ids(str(args.sample_file) if args.sample_file else None)]))

def _operon_definitions(args, run: Path) -> list[dict[str,object]]:
    if not args.manifest:
        if not args.genes: raise SystemExit("get-operon requires --genes or --manifest")
        samples=_samples(args); organism=resolve_organism(run,args.organism,samples)
        return [{"name":args.operon_name or "","genes":args.genes,"organism":organism,"samples":samples}]
    definitions={}
    for row in read_manifest_rows(str(args.manifest)):
        name=row.get("operon_name") or row.get("operon") or ""
        genes=[x for x in re_split(row.get("genes") or row.get("gene") or "") if x]
        if not genes: raise SystemExit("Operon manifest requires gene or genes values")
        key=(name,row.get("organism", "")); item=definitions.setdefault(key,{"name":name,"genes":[],"organism":row.get("organism", ""),"samples":[]})
        item["genes"].extend(genes)
        if row.get("isolate_id"): item["samples"].append(row["isolate_id"])
    result=[]
    for item in definitions.values():
        item["genes"]=list(dict.fromkeys(item["genes"])); item["samples"]=list(dict.fromkeys(item["samples"])); item["organism"]=resolve_organism(run,item["organism"] or None,item["samples"]); result.append(item)
    return result

def re_split(value: str) -> list[str]:
    import re
    return [x.strip() for x in re.split(r"[,;\s]+",value) if x.strip()]

def _differential_cohorts(args, run: Path, organism: str) -> tuple[list[str],list[str]]:
    a=[*args.cohort_a,*read_ids(str(args.cohort_a_file) if args.cohort_a_file else None),*_samples(args)]; b=[*args.cohort_b,*read_ids(str(args.cohort_b_file) if args.cohort_b_file else None)]
    if args.manifest:
        rows=read_manifest_rows(str(args.manifest)); labels=[]
        if any(not row.get("isolate_id") for row in rows): raise SystemExit("Differential manifest requires isolate_id for every row")
        for row in rows:
            if row.get(args.group_column) and row[args.group_column] not in labels: labels.append(row[args.group_column])
        if len(labels)!=2 and not (args.group_a_label and args.group_b_label): raise SystemExit(f"Differential manifest must contain exactly two {args.group_column} values, or provide --group-a-label/--group-b-label")
        la=args.group_a_label or labels[0]; lb=args.group_b_label or labels[1]
        a.extend(r["isolate_id"] for r in rows if r.get(args.group_column)==la); b.extend(r["isolate_id"] for r in rows if r.get(args.group_column)==lb)
        row_organisms={r.get("organism","") for r in rows if r.get("organism")}
        if row_organisms and row_organisms!={organism}: raise SystemExit("Differential manifest organism does not match the selected run organism")
    a=list(dict.fromkeys(a)); b=list(dict.fromkeys(b))
    if b and not a: raise SystemExit("A cohort B was provided without cohort A")
    return a,b

def _resolve_result(path: Path, filename: str) -> str:
    path=path.expanduser().resolve(); result=path/filename if path.is_dir() else path
    if not result.is_file(): raise SystemExit(f"Prior utility result not found: {result}")
    return str(result)

def _request(args, run: Path) -> dict[str,object]:
    samples=_samples(args); kind=args.utility.replace("-","_")
    if kind=="get_samples":
        organism=resolve_organism(run,args.organism,samples); validate_matrix_selection(run,organism,args.genes,samples)
        return {"utility":kind,"organism":organism,"genes":args.genes,"samples":samples,"status":args.status,"match":args.match}
    if kind=="get_differential_genes":
        manifest_rows=read_manifest_rows(str(args.manifest)) if args.manifest else []
        manifest_samples=[r["isolate_id"] for r in manifest_rows if r.get("isolate_id")]; manifest_organisms={r.get("organism","") for r in manifest_rows if r.get("organism")}
        if len(manifest_organisms)>1: raise SystemExit("Differential manifest spans multiple organisms")
        organism=resolve_organism(run,args.organism or (next(iter(manifest_organisms)) if manifest_organisms else None),samples+manifest_samples); a,b=_differential_cohorts(args,run,organism); validate_matrix_selection(run,organism,[],a+b)
        return {"utility":kind,"organism":organism,"group_a":a,"group_b":b,"max_q_value":args.max_q_value,"min_prevalence_difference":args.min_prevalence_difference,"top":args.top}
    if kind=="get_operon":
        definitions=_operon_definitions(args,run)
        for d in definitions:
            validate_matrix_selection(run,str(d["organism"]),list(d["genes"]),list(d["samples"]))
        return {"utility":kind,"operons":definitions}
    if kind=="get_variants":
        genes=list(args.genes or []); organism=args.organism
        if args.operon:
            calls=read_tsv(Path(_resolve_result(args.operon,"operon_calls.tsv"))); genes=list(dict.fromkeys(k for r in calls for k,v in r.items() if k not in {"operon_name","operon_id","organism","isolate_id"} and v!="")); samples=list(dict.fromkeys(r["isolate_id"] for r in calls)); organisms={r["organism"] for r in calls};
            if len(organisms)!=1: raise SystemExit("Selected operon result spans multiple organisms")
            organism=next(iter(organisms))
        if not genes: raise SystemExit("get-variants requires --genes or --operon")
        organism=resolve_organism(run,organism,samples); validate_matrix_selection(run,organism,genes,samples)
        if not 0<=args.min_similarity<=100: raise SystemExit("--min-similarity must be between 0 and 100")
        if not 0<=args.flanking_genes<=10: raise SystemExit("--flanking-genes must be between 0 and 10")
        return {"utility":kind,"organism":organism,"genes":genes,"samples":samples,"min_similarity":args.min_similarity,"flanking_genes":args.flanking_genes,"analyze_flanks":args.analyze_flanks}
    if kind=="diagnose_call":
        organism=resolve_organism(run,args.organism,samples); validate_matrix_selection(run,organism,args.genes,samples)
        if args.max_samples<1: raise SystemExit("--max-samples must be at least 1")
        return {"utility":kind,"organism":organism,"genes":args.genes,"samples":samples,"max_samples":args.max_samples,"replay_assembly":not args.skip_assembly_replay}
    if not args.genes and not args.operon and not args.variants: raise SystemExit("itol requires --genes, --operon, and/or --variants")
    inferred_organisms=set(); inferred_samples=[]
    for path,filename in ((args.operon,"operon_calls.tsv"),(args.variants,"gene_variants.tsv")):
        if not path: continue
        rows=read_tsv(Path(_resolve_result(path,filename))); inferred_organisms.update(r.get("organism","") for r in rows if r.get("organism")); inferred_samples.extend(r["isolate_id"] for r in rows if r.get("isolate_id"))
    if len(inferred_organisms)>1: raise SystemExit("Prior utility results span multiple organisms")
    if not samples and not args.organism: samples=list(dict.fromkeys(inferred_samples))
    organism=resolve_organism(run,args.organism or (next(iter(inferred_organisms)) if inferred_organisms else None),samples); validate_matrix_selection(run,organism,args.genes,samples)
    custom={}
    if args.custom_colors:
        custom={r["label"]:r["color"] for r in read_manifest_rows(str(args.custom_colors))}
    if args.color_scheme=="custom" and not custom: raise SystemExit("--color-scheme custom requires --custom-colors with label and color columns")
    return {"utility":"itol","organism":organism,"genes":args.genes,"samples":samples,"operon":_resolve_result(args.operon,"operon_calls.tsv") if args.operon else "","variants":_resolve_result(args.variants,"gene_variants.tsv") if args.variants else "","color_scheme":args.color_scheme,"custom_colors":custom}

def utils_command(args) -> int:
    print(welcome("Welcome to CleanGene Utils"))
    run=locate_run(args); print(f"Located run: {run}")
    with spinner("Getting ready to submit"):
        request=_request(args,run); cfg={**DEFAULTS,**json.loads((run/"provenance"/"resolved_config.json").read_text())}; stamp=datetime.now().strftime("%y%m%d_%H%M%S"); analysis_id=safe_name(args.analysis_name or f"{stamp}_{request['utility']}"); out=run/"results"/"utils"/analysis_id; logs=run/"logs"/"slurm"/"utils"; out.mkdir(parents=True,exist_ok=False); logs.mkdir(parents=True,exist_ok=True)
        resource="DIAGNOSTIC" if request["utility"]=="diagnose_call" else "VARIANT" if request["utility"]=="get_variants" else ""
        cpus=cfg.get(f"UTILS_{resource}_CPUS" if resource else "UTILS_CPUS","8"); mem=cfg.get(f"UTILS_{resource}_MEM" if resource else "UTILS_MEM","32G"); limit=cfg.get(f"UTILS_{resource}_TIME" if resource else "UTILS_TIME","12:00:00")
        request.update({"run_dir":str(run),"output_dir":str(out),"analysis_id":analysis_id,"submitted":datetime.now().isoformat(),"cpus":int(cpus),"min_mapq":int(cfg.get("READ_VALIDATION_MIN_MAPQ","20")),"min_depth":float(cfg.get("READ_VALIDATION_MIN_MEAN_DEPTH","5")),"basequal":int(cfg.get("BASEQUAL","30"))}); request_path=out/"request.json"; atomic_json(request_path,request)
        import shlex
        wrap=f"{shlex.quote(sys.executable)} -m cleangene _utils_worker --request {shlex.quote(str(request_path))}"
        command=sbatch_cmd(name=f"cg-util-{safe_name(str(request['utility']))[:20]}",wrap=wrap,cpus=cpus,mem=mem,time=limit,account=cfg.get("SLURM_ACCOUNT",""),partition=cfg.get("SLURM_PARTITION",""),log=logs/f"{analysis_id}.%j.log")
        job_id=submit(command,args.dry_run); request["slurm_job_id"]=job_id; atomic_json(request_path,request)
    print(submitted(f"Analysis submitted. Please find logs in {logs}"))
    return 0
