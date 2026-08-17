from __future__ import annotations
import csv, json, os, shutil, subprocess
from pathlib import Path
from .config import truthy
from .evidence import validate_isolate
from .fasta import assembly_metrics
from .manifest import load_manifest, groups
from .pangenome import normalize_panaroo, recover_sequences, select_rows, write_binary
from .util import command_exists, load_json, read_tsv, run, safe_name, touch_done, write_tsv

def context(run_dir: Path):
    cfg=load_json(run_dir/"provenance"/"resolved_config.json"); rows=read_tsv(run_dir/"provenance"/"manifest.tsv"); return cfg, rows

def task_row(run_dir: Path, kind: str, index: int) -> dict[str,str]:
    rows=read_tsv(run_dir/"state"/f"{kind}_tasks.tsv")
    if index<0 or index>=len(rows): raise SystemExit(f"Task index {index} outside {kind} task list")
    return rows[index]

def parse_kraken_report(path: Path, expected: str) -> tuple[str,float,float]:
    top=("",-1.0); contamination=0.0; unclassified=0.0; expected_norm=" ".join(expected.lower().split())
    if not path.is_file(): return "",0.0,0.0
    for line in path.read_text(errors="replace").splitlines():
        f=line.split("\t")
        if len(f)<6: continue
        pct=float(f[0]); rank=f[3].strip(); name=" ".join(f[5].strip().split()); norm=name.lower()
        if rank=="S":
            if pct>top[1]: top=(name,pct)
            if expected_norm and norm!=expected_norm: contamination += pct
        if rank=="U": unclassified=max(unclassified,pct)
    return top[0], contamination, max(0.0,100.0-unclassified)

def prepare_read_inputs(row: dict[str,str], out: Path, logs: Path, cfg: dict[str,str]) -> tuple[str,str,str,int]:
    r1=row.get("R1","").strip(); r2=row.get("R2","").strip(); mode=cfg.get("READ_TRIMMING_MODE","auto").strip().lower()
    method="manifest_fastq"; trimmed=0
    if row.get("raw_bam","").strip():
        reads=out/"reads"; reads.mkdir(exist_ok=True)
        r1=str(reads/"raw_R1.fastq.gz"); r2=str(reads/"raw_R2.fastq.gz")
        collated=reads/"raw_name_collated.bam"
        other=reads/"raw_unpaired.fastq.gz"
        threads=str(max(1,int(cfg.get("CPUS","4"))-1))
        run(["samtools","collate","-@",threads,"-o",str(collated),row["raw_bam"]],stdout=logs/"samtools-collate.stdout",stderr=logs/"samtools-collate.stderr")
        run(["samtools","fastq","-@",threads,"-1",r1,"-2",r2,"-0",str(other),"-s",str(other),"-n",str(collated)],
            stdout=logs/"samtools-fastq.stdout",stderr=logs/"samtools-fastq.stderr")
        method="raw_bam_samtools_fastq"
    if mode not in {"off","auto","always"}: raise SystemExit("READ_TRIMMING_MODE must be off, auto, or always")
    if mode in {"auto","always"}:
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
    return r1,r2,method,trimmed

def preprocess(run_dir: Path, index: int) -> None:
    cfg, _=context(run_dir); row=task_row(run_dir,"isolate",index); iso=row["isolate_id"]; group=row["group_id"]; safe=safe_name(iso)
    root=run_dir/"results"/"groups"/safe_name(group); out=root/"01_isolates"/safe; done=run_dir/"state"/"preprocess"/f"{safe}.done.json"
    if done.is_file():
        status=load_json(done)
        if status.get("excluded") or (out/"annotation"/f"{safe}.gff").is_file(): return
    out.mkdir(parents=True,exist_ok=True); logs=out/"logs"; logs.mkdir(exist_ok=True)
    r1,r2,read_method,adapter_trimmed=prepare_read_inputs(row,out,logs,cfg)
    expected=row.get("organism","").strip(); taxonomy=cfg.get("TAXONOMY_MODE","auto"); excluded=False; top=""; contam=0.0
    if taxonomy!="off":
        db=cfg.get("KRAKEN2_DB","")
        if not db: raise SystemExit("KRAKEN2_DB is required when TAXONOMY_MODE is not off")
        report=out/"kraken2.report.tsv"; output=out/"kraken2.output.tsv"
        run(["kraken2","--db",db,"--paired","--report",str(report),"--output",str(output),r1,r2],stdout=logs/"kraken2.stdout",stderr=logs/"kraken2.stderr")
        top,contam,_=parse_kraken_report(report,expected)
        if expected and contam>float(cfg["MAJOR_CONTAMINATION_THRESHOLD"]): excluded=True
    assembly=row.get("assembly","").strip()
    if excluded:
        metrics=assembly_metrics(Path(assembly)) if assembly else {"assembly_length":"","contigs":"","n50":"","l50":"","ambiguous_bases":"","gc_fraction":""}
        fields=["isolate_id","group_id","excluded","reason","top_species","contamination_pct","R1","R2","raw_bam","read_preprocessing","adapter_trimmed","assembly","assembly_length","contigs","n50","l50","ambiguous_bases","gc_fraction","gff"]
        data={"isolate_id":iso,"group_id":group,"excluded":1,"reason":"major_contamination","top_species":top,"contamination_pct":contam,"R1":r1,"R2":r2,"raw_bam":row.get("raw_bam",""),"read_preprocessing":read_method,"adapter_trimmed":adapter_trimmed,"assembly":assembly,"gff":"",**metrics}
        write_tsv(out/"qc.tsv",fields,[data]); touch_done(done,{"excluded":True,"reason":"major_contamination"}); return
    if not assembly:
        shov=out/"assembly"; shov.mkdir(exist_ok=True); assembly=str(shov/"contigs.fa")
        if not Path(assembly).is_file(): run(["shovill","--R1",r1,"--R2",r2,"--outdir",str(shov),"--cpus",cfg.get("CPUS","4"),"--force"],stdout=logs/"shovill.stdout",stderr=logs/"shovill.stderr")
    ann=out/"annotation"; gff=ann/f"{safe}.gff"
    if not gff.is_file():
        if ann.exists(): shutil.rmtree(ann)
        run(["prokka","--outdir",str(ann),"--prefix",safe,"--locustag",safe,"--cpus",cfg.get("CPUS","4"),"--force",assembly],stdout=logs/"prokka.stdout",stderr=logs/"prokka.stderr")
    metrics=assembly_metrics(Path(assembly)); fields=["isolate_id","group_id","excluded","reason","top_species","contamination_pct","R1","R2","raw_bam","read_preprocessing","adapter_trimmed","assembly","assembly_length","contigs","n50","l50","ambiguous_bases","gc_fraction","gff"]
    data={"isolate_id":iso,"group_id":group,"excluded":int(excluded),"reason":"major_contamination" if excluded else "","top_species":top,"contamination_pct":contam,"R1":r1,"R2":r2,"raw_bam":row.get("raw_bam",""),"read_preprocessing":read_method,"adapter_trimmed":adapter_trimmed,"assembly":assembly,"gff":str(gff),**metrics}
    write_tsv(out/"qc.tsv",fields,[data]); touch_done(done,{"excluded":excluded,"gff":str(gff)})

def retained_rows(run_dir: Path, group: str) -> list[dict[str,str]]:
    _, rows=context(run_dir); result=[]
    for row in rows:
        if row["group_id"]!=group: continue
        safe=safe_name(row["isolate_id"]); qc=run_dir/"results"/"groups"/safe_name(group)/"01_isolates"/safe/"qc.tsv"
        if not qc.is_file(): raise SystemExit(f"Missing isolate QC: {qc}")
        q=read_tsv(qc)[0]
        if q["excluded"] in {"1","true","True"}: continue
        row=dict(row); row["assembly"]=q["assembly"]; row["gff"]=q["gff"]; row["R1"]=q.get("R1",row.get("R1","")); row["R2"]=q.get("R2",row.get("R2","")); result.append(row)
    return result

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
    if done.is_file():
        status=load_json(done)
        if status.get("status")=="skipped" or status.get("external_pangenome_dir") or (out/"gene_presence_absence.csv").is_file(): return
    retained=retained_rows(run_dir,group)
    if len(retained)<2: touch_done(done,{"status":"skipped","reason":"fewer_than_two_retained"}); return
    external=manifest_pangenome_dir(rows,group)
    logs=root/"logs"; logs.mkdir(exist_ok=True)
    if external:
        isolates=[r["isolate_id"] for r in retained]; rows=normalize_panaroo(external/"gene_presence_absence.csv",isolates); calls=root/"02_pangenome"/"initial_calls"; calls.mkdir(parents=True,exist_ok=True); write_binary(calls/"gene_presence_absence.binary.tsv",rows,isolates)
        touch_done(done,{"n_isolates":len(isolates),"n_genes":len(rows),"external_pangenome_dir":str(external)})
        return
    out.mkdir(parents=True,exist_ok=True)
    gffs=[r["gff"] for r in retained]
    run(["panaroo","-i",*gffs,"-o",str(out),"--clean-mode",cfg["PANAROO_CLEAN_MODE"],"-t",cfg.get("PANAROO_CPUS",cfg.get("CPUS","4"))],stdout=logs/"panaroo.stdout",stderr=logs/"panaroo.stderr")
    isolates=[r["isolate_id"] for r in retained]; rows=normalize_panaroo(out/"gene_presence_absence.csv",isolates); calls=root/"02_pangenome"/"initial_calls"; calls.mkdir(parents=True,exist_ok=True); write_binary(calls/"gene_presence_absence.binary.tsv",rows,isolates)
    touch_done(done,{"n_isolates":len(isolates),"n_genes":len(rows)})

def prepare_validation(run_dir: Path, index: int) -> None:
    cfg,_=context(run_dir); group=task_row(run_dir,"group",index)["group_id"]; root=run_dir/"results"/"groups"/safe_name(group); done=run_dir/"state"/"prepare_validation"/f"{safe_name(group)}.done.json"
    if done.is_file(): return
    retained=retained_rows(run_dir,group)
    if len(retained)<2: touch_done(done,{"status":"skipped"}); return
    isolates=[r["isolate_id"] for r in retained]; initial=root/"02_pangenome"/"initial_calls"/"gene_presence_absence.binary.tsv"; panaroo_dir=prepared_pangenome_dir(run_dir,group,root)
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
        ev.mkdir(parents=True,exist_ok=True); write_tsv(ev/"metrics.tsv",["reference_id","Gene","validated_call","validation_state","breadth","percent_coverage","mean_depth","identity","percent_identity","identical_positions","aligned_positions","mapped_reads","mean_mapping_quality"],[]); touch_done(done,{"status":"no_selected_genes"}); return
    rr=retained[iso]; validate_isolate(ref,key,rr["R1"],rr["R2"],ev,int(cfg["VALIDATION_CPUS"]),float(cfg["READ_VALIDATION_MIN_BREADTH"]),float(cfg["READ_VALIDATION_MIN_MEAN_DEPTH"]),float(cfg["READ_VALIDATION_MIN_IDENTITY"]),int(cfg["READ_VALIDATION_MIN_MAPQ"]),int(cfg["BASEQUAL"])); touch_done(done)

def reduce_group(run_dir: Path, index: int) -> None:
    cfg,_=context(run_dir); group=task_row(run_dir,"group",index)["group_id"]; root=run_dir/"results"/"groups"/safe_name(group); done=run_dir/"state"/"reduce"/f"{safe_name(group)}.done.json"
    if done.is_file(): return
    retained=retained_rows(run_dir,group)
    if len(retained)<2: touch_done(done,{"status":"skipped"}); return
    isolates=[r["isolate_id"] for r in retained]; initial_path=root/"02_pangenome"/"initial_calls"/"gene_presence_absence.binary.tsv"; initial=read_tsv(initial_path); by_gene={r["Gene"]:{i:int(r[i]) for i in isolates} for r in initial}; out=root/"03_read_validation"; metrics=[]
    for iso in isolates:
        p=out/"evidence"/safe_name(iso)/"metrics.tsv"
        for r in read_tsv(p) if p.is_file() else []:
            r={**r,"isolate_id":iso,"initial_call":by_gene[r["Gene"]][iso]}
            r.setdefault("percent_coverage", str(float(r.get("breadth") or 0) * 100.0))
            r.setdefault("percent_identity", str(float(r.get("identity") or 0) * 100.0))
            metrics.append(r)
    validated={g:dict(v) for g,v in by_gene.items()}
    for r in metrics: validated[r["Gene"]][r["isolate_id"]]=int(r["validated_call"])
    fields=["Gene",*isolates]; write_tsv(out/"validated_gene_presence_absence.binary.tsv",fields,([g,*[validated[g][i] for i in isolates]] for g in by_gene))
    metric_fields=["Gene","isolate_id","initial_call","validated_call","validation_state","breadth","percent_coverage","mean_depth","identity","percent_identity","identical_positions","aligned_positions","mapped_reads","mean_mapping_quality"]
    write_tsv(out/"read_validation_metrics.tsv",metric_fields,metrics)
    evidence_rows=[]
    metric_index={(r["Gene"],r["isolate_id"]):r for r in metrics}
    for gene in by_gene:
        for iso in isolates:
            m=metric_index.get((gene,iso)); initial_call=by_gene[gene][iso]; final=validated[gene][iso]
            evidence_rows.append([group,gene,iso,initial_call,final,m["validation_state"] if m else "not_tested_carried_forward",m.get("breadth","") if m else "",m.get("percent_coverage","") if m else "",m.get("mean_depth","") if m else "",m.get("identity","") if m else "",m.get("percent_identity","") if m else "",m.get("mapped_reads","") if m else ""])
    write_tsv(out/"gene_call_evidence.long.tsv",["group_id","Gene","isolate_id","initial_call","validated_call","validation_state","breadth","percent_coverage","mean_depth","identity","percent_identity","mapped_reads"],evidence_rows)
    changes=[]
    for g in by_gene:
        a=[by_gene[g][i] for i in isolates]; b=[validated[g][i] for i in isolates]; changes.append([g,sum(a),sum(b),sum(x==0 and y==1 for x,y in zip(a,b)),sum(x==1 and y==0 for x,y in zip(a,b)),sum(x==y for x,y in zip(a,b))])
    write_tsv(out/"tested_genes.tsv",["Gene","n_initial_present","n_validated_present","n_added","n_removed","n_unchanged"],changes)
    touch_done(done,{"n_isolates":len(isolates),"n_genes":len(by_gene)})

def summarize(run_dir: Path) -> None:
    cfg,rows=context(run_dir); iso_rows=[]; group_rows=[]
    for group in groups(rows):
        root=run_dir/"results"/"groups"/safe_name(group); retained=retained_rows(run_dir,group); val=root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv"
        n_genes=max(0,len(read_tsv(val))) if val.is_file() else 0; group_rows.append([group,len([r for r in rows if r["group_id"]==group]),len(retained),n_genes])
        for row in rows:
            if row["group_id"]!=group: continue
            q=read_tsv(root/"01_isolates"/safe_name(row["isolate_id"])/"qc.tsv")[0]; iso_rows.append(q)
    cohort=run_dir/"results"/"cohort"; write_tsv(cohort/"isolate_qc.tsv",["isolate_id","group_id","excluded","reason","top_species","contamination_pct","R1","R2","raw_bam","read_preprocessing","adapter_trimmed","assembly","assembly_length","contigs","n50","l50","ambiguous_bases","gc_fraction","gff"],iso_rows); write_tsv(cohort/"group_summary.tsv",["group_id","input_isolates","retained_isolates","validated_gene_clusters"],group_rows); touch_done(run_dir/"state"/"summary.done.json",{"groups":len(group_rows)})

def dispatch(stage: str, run_dir: Path, index: int | None) -> None:
    if stage=="preprocess": preprocess(run_dir,int(index))
    elif stage=="panaroo": panaroo(run_dir,int(index))
    elif stage=="prepare_validation": prepare_validation(run_dir,int(index))
    elif stage=="validate": validate(run_dir,int(index))
    elif stage=="reduce": reduce_group(run_dir,int(index))
    elif stage=="summary": summarize(run_dir)
    else: raise SystemExit(f"Unknown worker stage: {stage}")
