from __future__ import annotations
import csv, gzip, shutil, subprocess
from pathlib import Path
from .fasta import write_fasta
from .util import read_tsv, run, safe_name, write_tsv

def _matrix(path: Path) -> tuple[list[str],dict[str,dict[str,int]]]:
    rows=read_tsv(path); isolates=list(rows[0])[1:] if rows else []
    return isolates,{row["Gene"]:{iso:int(row[iso]) for iso in isolates} for row in rows}

def _norm_read(name: str) -> str:
    return name.split()[0].removesuffix("/1").removesuffix("/2")

def _bam_records(path: Path, source: str, ref_to_gene: dict[str,str]) -> list[dict[str,object]]:
    if not path.is_file(): return []
    result=subprocess.run(["samtools","view",str(path)],check=True,capture_output=True,text=True); rows=[]
    for line in result.stdout.splitlines():
        f=line.split("\t")
        if len(f)<11 or f[2] not in ref_to_gene: continue
        rows.append({"isolate_id":"","Gene":ref_to_gene[f[2]],"read_id":_norm_read(f[0]),"source":source,"reference_id":f[2],"flag":f[1],"position":f[3],"mapq":f[4],"cigar":f[5],"sequence":f[9]})
    return rows

def _names_by_gene(rows: list[dict[str,object]]) -> dict[str,set[str]]:
    result={}
    for row in rows: result.setdefault(str(row["Gene"]),set()).add(str(row["read_id"]))
    return result

def _open_fastq(path: Path):
    return gzip.open(path,"rt",errors="replace") if path.suffix==".gz" else path.open(errors="replace")

def _fastq_members(paths: list[Path], wanted: set[str]) -> set[str]:
    found=set()
    for path in paths:
        if not path.is_file(): continue
        with _open_fastq(path) as handle:
            while True:
                name=handle.readline()
                if not name: break
                handle.readline(); handle.readline(); handle.readline()
                key=_norm_read(name[1:].strip())
                if key in wanted: found.add(key)
        if found==wanted: break
    return found

def _map_reads(reference: Path, r1: str, r2: str, out: Path, threads: int, min_mapq: int) -> Path:
    from .evidence import map_reads
    out.parent.mkdir(parents=True,exist_ok=True); map_reads(reference,r1,r2,out,threads,min_mapq,out.with_suffix(".bwa.log")); return out

def _target_support(target: Path | None, label: str, r1: str, r2: str, wanted: set[str], out: Path, threads: int, min_mapq: int) -> tuple[set[str],str]:
    if not target or not target.is_file(): return set(),"unavailable"
    local=out/f"{label}.fasta"
    local.parent.mkdir(parents=True,exist_ok=True)
    if target.suffix==".gz":
        with gzip.open(target,"rb") as source, local.open("wb") as sink: shutil.copyfileobj(source,sink)
    else:
        shutil.copy2(target,local)
    run(["bwa","index",str(local)]); bam=_map_reads(local,r1,r2,out/f"{label}.bam",threads,min_mapq)
    result=subprocess.run(["samtools","view",str(bam)],check=True,capture_output=True,text=True)
    mapped={_norm_read(line.split("\t",1)[0]) for line in result.stdout.splitlines() if line}
    return wanted&mapped,str(target)

def _first(root: Path, names: tuple[str,...]) -> Path | None:
    for name in names:
        hits=list(root.rglob(name))
        if hits: return hits[0]
    return None

def _tool_versions(out: Path) -> None:
    rows=[]
    for tool,args in (("shovill",["--version"]),("kmc",["-h"]),("seqkit",["version"]),("spades.py",["--version"]),("fastp",["--version"]),("prokka",["--version"]),("panaroo",["--version"]),("bwa",[]),("samtools",["--version"])):
        try:
            result=subprocess.run([tool,*args],capture_output=True,text=True,timeout=30); text=(result.stdout or result.stderr).strip().splitlines(); rows.append([tool,text[0] if text else "available; version not reported"])
        except (OSError,subprocess.TimeoutExpired): rows.append([tool,"unavailable"])
    write_tsv(out/"tool_versions.tsv",["program","version"],rows)

def infer_failure_stage(initial: int,validated: int,raw: int,processed: int,sampled: int,graph: int,spades: int,replay_final: int,final_locus: bool,annotated: bool,replay: bool=True) -> tuple[str,str]:
    if initial==1 and validated==0: return "read_validation","BWA/SAMtools evidence rejected an initial Panaroo presence call"
    if initial==validated: return "no_call_change","The selected call was not changed by validation"
    if raw>processed: return "fastp","Gene-supporting raw reads were absent after read preprocessing"
    if replay and processed>sampled: return "shovill_seqkit_sampling","Gene-supporting reads were present after trimming but absent from Shovill depth sampling"
    if replay and sampled and graph==0: return "spades_graph_construction","Sampled supporting reads were available to SPAdes but did not map to its retained assembly graph"
    if replay and graph and spades==0: return "spades_contig_resolution","Supporting reads mapped to the SPAdes graph but not its contig output"
    if replay and spades and replay_final==0: return "shovill_post_assembly","SPAdes contig support was lost during Shovill correction or contig filtering"
    if final_locus and not annotated: return "prokka","The gene sequence was found on the final contigs but no overlapping Prokka feature was found"
    if final_locus and initial==0: return "panaroo","The gene sequence and annotation context reached the final assembly, but Panaroo made an absence call"
    return "unresolved","Available artifacts do not isolate one responsible program"

def _annotation_overlap(reference: Path, assembly: Path, gff: Path) -> dict[str,bool]:
    from .downstream import _gff_features, _locate
    locations=_locate(reference,assembly) if assembly.is_file() else {}; features=_gff_features(gff) if gff.is_file() else []
    result={}
    for ref,loc in locations.items():
        result[ref]=any(f["contig"]==loc["contig"] and int(f["start"])<=int(loc["end"]) and int(f["end"])>=int(loc["start"]) for f in features)
    return result

def diagnose_call(request: dict[str,object]) -> None:
    from .downstream import _gene_references
    from .workers import find_isolate_qc, retained_rows
    run_dir=Path(str(request["run_dir"])); out=Path(str(request["output_dir"])); organism=str(request["organism"]); genes=list(request["genes"]); threads=int(request.get("cpus",8)); min_mapq=int(request.get("min_mapq",20)); replay=bool(request.get("replay_assembly",True))
    root=run_dir/"results"/"groups"/safe_name(organism); _,initial=_matrix(root/"02_pangenome"/"initial_calls"/"gene_presence_absence.binary.tsv"); isolates,validated=_matrix(root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv")
    requested=list(request.get("samples",[])); discordant=[iso for iso in isolates if any(initial[g][iso]!=validated[g][iso] for g in genes)]
    selected=requested or discordant[:int(request.get("max_samples",20))]
    if not selected: raise SystemExit("No initial-versus-validated call changes matched the selected genes; provide --samples to diagnose unchanged calls")
    retained={row["isolate_id"]:row for row in retained_rows(run_dir,organism)}; manifest={row["isolate_id"]:row for row in read_tsv(run_dir/"provenance"/"manifest.tsv")}
    missing=[iso for iso in selected if iso not in retained]
    if missing: raise SystemExit("Diagnostic samples are not retained in this organism: "+", ".join(missing[:20]))
    records,key=_gene_references(run_dir,organism,genes); reference=out/"diagnostic_gene_references.fasta"; write_fasta(reference,records); run(["bwa","index",str(reference)]); ref_to_gene={str(row["reference_id"]):str(row["Gene"]) for row in key}; refs_by_gene={g:{r for r,x in ref_to_gene.items() if x==g} for g in genes}
    metric_rows=read_tsv(root/"03_read_validation"/"read_validation_metrics.tsv"); metrics={(r["isolate_id"],r["Gene"]):r for r in metric_rows}; summaries=[]; read_rows=[]; artifacts=[]; _tool_versions(out)
    for iso in selected:
        sample=out/"samples"/safe_name(iso); sample.mkdir(parents=True,exist_ok=True); row=retained[iso]; source=manifest[iso]; processed_r1=row["R1"]; processed_r2=row["R2"]
        validation_bam=root/"03_read_validation"/"evidence"/safe_name(iso)/"gene_reads.bam"; validation_rows=_bam_records(validation_bam,"validation_bwa",ref_to_gene); validation_names=_names_by_gene(validation_rows)
        raw_r1=source.get("R1",""); raw_r2=source.get("R2",""); qc_dir=find_isolate_qc(run_dir,source).parent
        if source.get("raw_bam") and (qc_dir/"reads"/"raw_R1.fastq.gz").is_file() and (qc_dir/"reads"/"raw_R2.fastq.gz").is_file():
            raw_r1=str(qc_dir/"reads"/"raw_R1.fastq.gz"); raw_r2=str(qc_dir/"reads"/"raw_R2.fastq.gz")
        raw_rows=[]
        if raw_r1 and raw_r2 and (raw_r1,raw_r2)!=(processed_r1,processed_r2):
            raw_rows=_bam_records(_map_reads(reference,raw_r1,raw_r2,sample/"raw_gene_reads.bam",threads,min_mapq),"raw_bwa",ref_to_gene)
        else: raw_rows=[{**r,"source":"raw_same_as_processed"} for r in validation_rows]
        processed_rows=_bam_records(_map_reads(reference,processed_r1,processed_r2,sample/"processed_gene_reads.bam",threads,min_mapq),"processed_replay_bwa",ref_to_gene)
        raw_names=_names_by_gene(raw_rows); processed_names=_names_by_gene(processed_rows); read_rows.extend({**r,"isolate_id":iso} for r in [*validation_rows,*raw_rows,*processed_rows])
        replay_dir=sample/"shovill_replay"; replay_status="skipped"; sampled_files=[]
        if replay:
            command=["shovill","--R1",processed_r1,"--R2",processed_r2,"--outdir",str(replay_dir),"--tmpdir",str(sample/"shovill_tmp"),"--cpus",str(threads),"--keepfiles","--force"]
            try: run(command,stdout=sample/"shovill_replay.stdout",stderr=sample/"shovill_replay.stderr"); replay_status="complete"
            except subprocess.CalledProcessError: replay_status="failed"
            sampled_files=list(sample.rglob("*.sub.fq.gz"))
        graph=_first(sample,("assembly_graph.fastg","assembly_graph.fasta")); spades=_first(sample,("spades.fasta","contigs.fasta")); replay_final=_first(replay_dir,("contigs.fa",)); original_assembly=Path(row.get("assembly","") or "missing"); gff=Path(row.get("gff","") or "missing"); annotation=_annotation_overlap(reference,original_assembly,gff)
        all_wanted={name for source_names in (raw_names,processed_names,validation_names) for values in source_names.values() for name in values}; sampled_members=_fastq_members(sampled_files,all_wanted) if sampled_files else (all_wanted if replay_status=="complete" else set())
        target_counts={}
        for label,target in (("spades_graph",graph),("spades_contigs",spades),("shovill_final",replay_final),("original_final",original_assembly)):
            target_counts[label]=_target_support(target,label,processed_r1,processed_r2,all_wanted,sample/"target_mappings",threads,min_mapq)
        for gene in genes:
            raw_gene=raw_names.get(gene,set()); processed_gene=processed_names.get(gene,set()); wanted=validation_names.get(gene,set()); sampled_gene=processed_gene&sampled_members; sampled=len(sampled_gene)
            graph_n=len(sampled_gene&target_counts["spades_graph"][0]); spades_n=len(sampled_gene&target_counts["spades_contigs"][0]); replay_n=len(sampled_gene&target_counts["shovill_final"][0]); original_n=len(processed_gene&target_counts["original_final"][0])
            gene_locations=set(refs_by_gene[gene])&set(annotation); final_locus=bool(gene_locations); annotated=any(annotation[r] for r in gene_locations)
            stage,reason=infer_failure_stage(initial[gene][iso],validated[gene][iso],len(raw_gene),len(processed_gene),sampled,graph_n,spades_n,replay_n,final_locus,annotated,replay_status=="complete")
            metric=metrics.get((iso,gene),{}); summaries.append([iso,organism,gene,initial[gene][iso],validated[gene][iso],metric.get("validation_state",""),metric.get("percent_coverage",""),metric.get("percent_identity",""),len(raw_gene),len(processed_gene),len(wanted),sampled,graph_n,spades_n,replay_n,original_n,int(final_locus),int(annotated),replay_status,stage,reason,"KMC estimates k-mer/genome statistics; Shovill's seqkit depth sampling is the direct read-removal step"])
        for label,target in (("validation_bam",validation_bam),("processed_R1",Path(processed_r1)),("processed_R2",Path(processed_r2)),("original_assembly",original_assembly),("prokka_gff",gff),("spades_graph",graph),("spades_contigs",spades),("shovill_replay_final",replay_final)):
            artifacts.append([iso,label,"" if target is None else str(target),int(bool(target and target.exists()))])
    fields=["isolate_id","organism","Gene","initial_call","validated_call","validation_state","percent_coverage","percent_identity","raw_support_reads","processed_support_reads","validation_support_reads","support_reads_retained_after_sampling","support_reads_mapping_spades_graph","support_reads_mapping_spades_contigs","support_reads_mapping_shovill_replay_final","support_reads_mapping_original_final","gene_sequence_on_original_final","overlapping_prokka_feature","replay_status","likely_failure_stage","interpretation","kmc_assessment"]
    write_tsv(out/"diagnostic_summary.tsv",fields,summaries); write_tsv(out/"supporting_reads.tsv",["isolate_id","Gene","read_id","source","reference_id","flag","position","mapq","cigar","sequence"],read_rows); write_tsv(out/"source_artifacts.tsv",["isolate_id","artifact","path","exists"],artifacts)
