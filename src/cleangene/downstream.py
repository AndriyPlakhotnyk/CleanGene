from __future__ import annotations
import csv, hashlib, math, subprocess
from collections import Counter
from pathlib import Path
from .fasta import read_fasta, write_fasta
from .pangenome import recover_sequences
from .util import load_json, read_tsv, run, safe_name, write_tsv

def available_organisms(run_dir: Path) -> list[str]:
    tasks=run_dir/"state"/"group_tasks.tsv"
    return [r["group_id"] for r in read_tsv(tasks)] if tasks.is_file() else []

def resolve_organism(run_dir: Path, organism: str | None, samples: list[str] | None = None) -> str:
    organisms=available_organisms(run_dir)
    if organism:
        if organism not in organisms:
            raise SystemExit(f"Organism not found in run: {organism}\nAvailable organisms: {', '.join(organisms) or '<none>'}")
        return organism
    if samples:
        by_iso={r["isolate_id"]:r["group_id"] for r in read_tsv(run_dir/"provenance"/"manifest.tsv")}
        matched={by_iso[x] for x in samples if x in by_iso}
        missing=[x for x in samples if x not in by_iso]
        if missing: raise SystemExit("Sample IDs not found in run: " + ", ".join(missing[:20]))
        if len(matched)==1: return next(iter(matched))
        if len(matched)>1: raise SystemExit("Samples span multiple organisms; provide --organism")
    raise SystemExit("Organism was not provided. Available organisms: " + (", ".join(organisms) or "<none>"))

def matrix_path(run_dir: Path, organism: str) -> Path:
    root=run_dir/"results"/"groups"/safe_name(organism); primary=root/"cleaned_pangenome.tsv"
    path=primary if primary.is_file() else root/"03_read_validation"/"validated_gene_presence_absence.binary.tsv"
    if not path.is_file(): raise SystemExit(f"Validated presence/absence matrix not found for {organism}: {path}")
    return path

def load_matrix(run_dir: Path, organism: str) -> tuple[list[str],dict[str,dict[str,int]]]:
    rows=read_tsv(matrix_path(run_dir,organism)); isolates=list(rows[0].keys())[1:] if rows else []
    return isolates,{r["Gene"]:{iso:int(r[iso]) for iso in isolates} for r in rows}

def select_genes(matrix: dict[str,dict[str,int]], genes: list[str]) -> list[str]:
    requested=list(dict.fromkeys(genes)); missing=[g for g in requested if g not in matrix]
    if missing: raise SystemExit("Genes not found in validated matrix: " + ", ".join(missing[:20]))
    return requested

def validate_matrix_selection(run_dir: Path, organism: str, genes: list[str], samples: list[str]) -> None:
    path=matrix_path(run_dir,organism)
    with path.open(newline="",errors="replace") as h:
        reader=csv.reader(h,delimiter="\t"); header=next(reader,[]); isolates=set(header[1:]); missing_samples=[x for x in samples if x not in isolates]
        if missing_samples: raise SystemExit("Samples not found in organism matrix: " + ", ".join(missing_samples[:20]))
        wanted=set(genes); found=set()
        for row in reader:
            if row and row[0] in wanted:
                found.add(row[0])
                if found==wanted: break
    missing_genes=[x for x in genes if x not in found]
    if missing_genes: raise SystemExit("Genes not found in validated matrix: " + ", ".join(missing_genes[:20]))

def read_ids(path: str | None) -> list[str]:
    if not path: return []
    values=[x.strip().split("\t")[0] for x in Path(path).read_text().splitlines() if x.strip() and not x.lstrip().startswith("#")]
    return values[1:] if values and values[0].lower() in {"isolate_id","sample_id"} else values

def read_manifest_rows(path: str) -> list[dict[str,str]]:
    p=Path(path)
    if not p.is_file(): raise SystemExit(f"Utility manifest not found: {p}")
    with p.open(newline="",encoding="utf-8-sig",errors="replace") as h:
        lines=[x for x in h if x.strip() and not x.lstrip().startswith("#")]
    if not lines: raise SystemExit(f"Utility manifest is empty: {p}")
    delimiter="," if lines[0].count(",")>lines[0].count("\t") else "\t"
    return [{(k or "").strip():(v or "").strip() for k,v in r.items()} for r in csv.DictReader(lines,delimiter=delimiter)]

def subset_samples(all_isolates: list[str], requested: list[str]) -> list[str]:
    if not requested: return all_isolates
    known=set(all_isolates); missing=[x for x in requested if x not in known]
    if missing: raise SystemExit("Samples not found in organism matrix: " + ", ".join(missing[:20]))
    selected=set(requested)
    return [x for x in all_isolates if x in selected]

def plot_binary_heatmap(path: Path, title: str, genes: list[str], isolates: list[str], values: list[list[int]]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    width=max(7,min(24,2+len(isolates)*0.12)); height=max(3,min(20,1.5+len(genes)*0.28))
    fig,ax=plt.subplots(figsize=(width,height)); ax.imshow(values,aspect="auto",interpolation="nearest",cmap=ListedColormap(["#f4f1ea","#263238"]))
    ax.set_title(title); ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes,fontsize=7)
    if len(isolates)<=100: ax.set_xticks(range(len(isolates))); ax.set_xticklabels(isolates,rotation=90,fontsize=6)
    else: ax.set_xticks([])
    ax.set_xlabel("isolates"); fig.tight_layout(); fig.savefig(path.with_suffix(".svg")); fig.savefig(path.with_suffix(".png"),dpi=180); plt.close(fig)

def get_samples(request: dict[str,object]) -> None:
    run_dir=Path(str(request["run_dir"])); out=Path(str(request["output_dir"])); organism=str(request["organism"]); genes=list(request["genes"])
    isolates,matrix=load_matrix(run_dir,organism); genes=select_genes(matrix,genes); requested=list(request.get("samples",[])); isolates=subset_samples(isolates,requested)
    fields=["isolate_id",*genes]; rows=[]; classified=[]
    for iso in isolates:
        calls=[matrix[g][iso] for g in genes]; rows.append([iso,*calls]); classified.append((iso,all(calls),any(calls)))
    write_tsv(out/"gene_presence_absence.tsv",fields,rows)
    mode=str(request.get("match","all")); status=str(request.get("status","both"))
    predicate=(lambda r:r[1]) if mode=="all" else (lambda r:r[2]); present=[[r[0]] for r in classified if predicate(r)]; absent=[[r[0]] for r in classified if not predicate(r)]
    if status in {"present","both"}: write_tsv(out/"samples_present.tsv",["isolate_id"],present)
    if status in {"absent","both"}: write_tsv(out/"samples_absent.tsv",["isolate_id"],absent)

def _fisher(a: int,b: int,c: int,d: int) -> tuple[float,float]:
    try:
        from scipy.stats import fisher_exact
        result=fisher_exact([[a,b],[c,d]],alternative="two-sided")
        return float(result.statistic),float(result.pvalue)
    except ImportError:
        n=a+b+c+d; row=a+b; col=a+c; lo=max(0,row-(n-col)); hi=min(row,col)
        def prob(x:int)->float: return math.comb(col,x)*math.comb(n-col,row-x)/math.comb(n,row)
        observed=prob(a); p=sum(prob(x) for x in range(lo,hi+1) if prob(x)<=observed+1e-15)
        return ((a*d)/(b*c) if b*c else math.inf if a*d else 0.0),min(1.0,p)

def _bh(pvalues: list[float]) -> list[float]:
    result=[1.0]*len(pvalues); previous=1.0
    for rank,i in reversed(list(enumerate(sorted(range(len(pvalues)),key=pvalues.__getitem__),1))):
        previous=min(previous,pvalues[i]*len(pvalues)/rank); result[i]=previous
    return result

def differential_genes(request: dict[str,object]) -> None:
    run_dir=Path(str(request["run_dir"])); out=Path(str(request["output_dir"])); organism=str(request["organism"])
    isolates,matrix=load_matrix(run_dir,organism); group_a=list(request.get("group_a",[])); group_b=list(request.get("group_b",[]))
    if not group_a:
        rows=[]
        for gene,calls in matrix.items():
            present=sum(calls.values()); prevalence=present/len(isolates) if isolates else 0
            rows.append([gene,present,len(isolates)-present,prevalence,1-abs(prevalence-0.5)*2])
        rows.sort(key=lambda x:(-x[-1],x[0])); write_tsv(out/"differential_genes.tsv",["Gene","present","absent","prevalence","variability_score"],rows)
        chosen=[r[0] for r in rows[:int(request.get("top",50))]]
        plot_binary_heatmap(out/"differential_genes_heatmap",f"{organism}: variable genes",chosen,isolates,[[matrix[g][i] for i in isolates] for g in chosen]); return
    group_a=subset_samples(isolates,group_a); group_b=subset_samples(isolates,group_b or [x for x in isolates if x not in set(group_a)])
    overlap=sorted(set(group_a)&set(group_b))
    if overlap: raise SystemExit("Differential cohorts overlap: " + ", ".join(overlap[:20]))
    if not group_a or not group_b: raise SystemExit("Differential analysis requires two non-empty cohorts")
    raw=[]
    for gene,calls in matrix.items():
        ap=sum(calls[i] for i in group_a); bp=sum(calls[i] for i in group_b); odds,p=_fisher(ap,len(group_a)-ap,bp,len(group_b)-bp)
        raw.append([gene,ap,len(group_a)-ap,bp,len(group_b)-bp,ap/len(group_a),bp/len(group_b),ap/len(group_a)-bp/len(group_b),odds,p])
    q=_bh([r[-1] for r in raw]); rows=[[*r,qv] for r,qv in zip(raw,q)]
    rows.sort(key=lambda r:(r[-1],-abs(r[7]),r[0])); fields=["Gene","group_a_present","group_a_absent","group_b_present","group_b_absent","group_a_prevalence","group_b_prevalence","prevalence_difference","odds_ratio","p_value","q_value"]
    write_tsv(out/"differential_genes.tsv",fields,rows)
    cutoff=float(request.get("max_q_value",0.05)); min_diff=float(request.get("min_prevalence_difference",0.1)); selected=[r[0] for r in rows if r[-1]<=cutoff and abs(r[7])>=min_diff][:int(request.get("top",50))]
    if not selected: selected=[r[0] for r in rows[:int(request.get("top",50))]]
    ordered=group_a+group_b; plot_binary_heatmap(out/"differential_genes_heatmap",f"{organism}: differential genes",selected,ordered,[[matrix[g][i] for i in ordered] for g in selected])
    write_tsv(out/"cohorts.tsv",["isolate_id","cohort"],([[x,"group_a"] for x in group_a]+[[x,"group_b"] for x in group_b]))

def get_operon(request: dict[str,object]) -> None:
    run_dir=Path(str(request["run_dir"])); out=Path(str(request["output_dir"])); definitions=list(request["operons"]); all_calls=[]; summaries=[]
    for number,definition in enumerate(definitions,1):
        organism=str(definition["organism"]); genes=list(definition["genes"]); name=str(definition.get("name","")).strip(); requested=list(definition.get("samples",[]))
        isolates,matrix=load_matrix(run_dir,organism); genes=select_genes(matrix,genes); isolates=subset_samples(isolates,requested)
        patterns=sorted({tuple(matrix[g][iso] for g in genes) for iso in isolates},reverse=True)
        ids={p:(f"{safe_name(name)}_V{i:03d}" if name else str(i)) for i,p in enumerate(patterns,1)}
        for iso in isolates:
            pattern=tuple(matrix[g][iso] for g in genes); all_calls.append({"operon_name":name or f"operon_{number}","operon_id":ids[pattern],"organism":organism,"isolate_id":iso,**{g:matrix[g][iso] for g in genes}})
        counts=Counter(ids[tuple(matrix[g][iso] for g in genes)] for iso in isolates)
        for pattern,oid in ids.items(): summaries.append({"operon_name":name or f"operon_{number}","operon_id":oid,"organism":organism,"pattern":"".join(map(str,pattern)),"genes":";".join(genes),"n_isolates":counts[oid]})
        plot_binary_heatmap(out/f"{safe_name(name or f'operon_{number}')}_heatmap",f"{name or f'Operon {number}'}: {organism}",genes,isolates,[[matrix[g][i] for i in isolates] for g in genes])
    gene_fields=sorted({g for d in definitions for g in d["genes"]})
    write_tsv(out/"operon_calls.tsv",["operon_name","operon_id","organism","isolate_id",*gene_fields],all_calls)
    write_tsv(out/"operon_summary.tsv",["operon_name","operon_id","organism","pattern","genes","n_isolates"],summaries)

def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]

def _gff_features(path: Path) -> list[dict[str,object]]:
    features=[]
    if not path.is_file(): return features
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("##FASTA"): break
        if not line or line.startswith("#"): continue
        f=line.split("\t")
        if len(f)<9 or f[2]!="CDS": continue
        attrs={k:v for k,v in (x.split("=",1) for x in f[8].split(";") if "=" in x)}
        from urllib.parse import unquote
        name=unquote(attrs.get("gene") or attrs.get("product") or attrs.get("locus_tag") or attrs.get("ID") or "unknown")
        features.append({"contig":f[0],"start":int(f[3]),"end":int(f[4]),"strand":f[6],"name":name,"locus_tag":attrs.get("locus_tag",attrs.get("ID",""))})
    return features

def _locate(reference: Path, assembly: Path) -> dict[str,dict[str,object]]:
    if not assembly.is_file(): return {}
    p=subprocess.run(["minimap2","-x","asm5","--secondary=no","-c",str(assembly),str(reference)],check=True,capture_output=True,text=True)
    located={}
    for line in p.stdout.splitlines():
        f=line.split("\t")
        if len(f)<12: continue
        candidate={"contig":f[5],"start":int(f[7])+1,"end":int(f[8]),"strand":f[4],"matches":int(f[9]),"block":int(f[10]),"mapq":int(f[11])}
        if f[0] not in located or candidate["matches"]>located[f[0]]["matches"]: located[f[0]]=candidate
    return located

def _neighbor_features(features: list[dict[str,object]], location: dict[str,object], count: int) -> list[tuple[int,dict[str,object]]]:
    same=sorted((x for x in features if x["contig"]==location["contig"]),key=lambda x:(x["start"],x["end"]))
    if not same: return []
    center=min(range(len(same)),key=lambda i:abs((int(same[i]["start"])+int(same[i]["end"]))/2-(int(location["start"])+int(location["end"]))/2))
    result=[]
    for offset in range(-count,count+1):
        if offset and 0<=center+offset<len(same): result.append((offset,same[center+offset]))
    return result

def _feature_sequence(assembly: Path, feature: dict[str,object]) -> str:
    seq=read_fasta(assembly).get(str(feature["contig"]),"")[int(feature["start"])-1:int(feature["end"])]
    return _revcomp(seq) if feature["strand"]=="-" else seq

def _variant_counts(vcf: Path) -> dict[str,dict[str,int]]:
    result: dict[str,dict[str,int]]={}
    if not vcf.is_file(): return result
    p=subprocess.run(["bcftools","query","-f","%CHROM\\t%POS\\t%REF\\t%ALT\\n",str(vcf)],check=True,capture_output=True,text=True)
    for line in p.stdout.splitlines():
        chrom,_,ref,alt=line.split("\t")[:4]; alt=alt.split(",")[0]; counts=result.setdefault(chrom,{"snps":0,"mnps":0,"insertions":0,"deletions":0,"inserted_bases":0,"deleted_bases":0})
        if len(ref)==len(alt)==1: counts["snps"]+=1
        elif len(ref)==len(alt): counts["mnps"]+=sum(a!=b for a,b in zip(ref,alt))
        elif len(alt)>len(ref): counts["insertions"]+=1; counts["inserted_bases"]+=len(alt)-len(ref)
        else: counts["deletions"]+=1; counts["deleted_bases"]+=len(ref)-len(alt)
    return result

def _plot_variant_alignment(out: Path, records: list[tuple[str,str]], counts: dict[str,int]) -> None:
    raw=out/"unique_variants.fasta"; aligned=out/"unique_variants.aligned.fasta"; write_fasta(raw,records)
    if not records:
        aligned.write_text("")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,3)); ax.text(0.5,0.5,"No variants met the minimum similarity",ha="center",va="center"); ax.axis("off"); fig.tight_layout(); fig.savefig(out/"unique_variant_alignment.pdf"); plt.close(fig); return
    try:
        with aligned.open("w") as h: subprocess.run(["mafft","--auto","--quiet",str(raw)],check=True,stdout=h)
        sequences=read_fasta(aligned)
    except (FileNotFoundError,subprocess.CalledProcessError):
        width=max(len(x[1]) for x in records); sequences={k:v.ljust(width,"-") for k,v in records}
        write_fasta(aligned,list(sequences.items()))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap
    names=list(sequences); palette=ListedColormap(["#4daf4a","#377eb8","#ff7f00","#e41a1c","#f2f2f2","#222222"]); code={"A":0,"C":1,"G":2,"T":3,"-":4,"N":5}
    with PdfPages(out/"unique_variant_alignment.pdf") as pdf:
        for start in range(0,len(names),80):
            page=names[start:start+80]; data=np.array([[code.get(x,5) for x in sequences[name]] for name in page],dtype=int)
            fig,ax=plt.subplots(figsize=(16,max(3,len(page)*0.24+1))); ax.imshow(data,aspect="auto",interpolation="nearest",cmap=palette,vmin=0,vmax=5,rasterized=True)
            ax.set_yticks(range(len(page))); ax.set_yticklabels([f"{x} (n={counts[x]})" for x in page],fontsize=7); ax.set_xlabel("aligned nucleotide position"); ax.set_title("Unique gene variants")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

def _gene_references(run_dir: Path, organism: str, genes: list[str]) -> tuple[list[tuple[str,str]],list[dict[str,str]]]:
    isolates,matrix=load_matrix(run_dir,organism); genes=select_genes(matrix,genes)
    root=run_dir/"results"/"groups"/safe_name(organism)
    from .workers import prepared_pangenome_dir
    panaroo=prepared_pangenome_dir(run_dir,organism,root)
    records,sources=recover_sequences([{"Gene":g} for g in genes],panaroo)
    return records,[{"reference_id":str(r[0]),"Gene":str(r[1]),"feature_type":"target","parent_gene":"","flank_offset":"","annotation":str(r[1])} for r in sources]

def get_variants(request: dict[str,object]) -> None:
    from .evidence import align_identity, consensus, coverage, fixed_coordinate_identity, map_reads
    from .workers import retained_rows
    run_dir=Path(str(request["run_dir"])); out=Path(str(request["output_dir"])); organism=str(request["organism"]); genes=list(request["genes"]); threads=int(request.get("cpus",4)); flank_count=int(request.get("flanking_genes",0)); analyze_flanks=bool(request.get("analyze_flanks",False)); min_similarity=float(request.get("min_similarity",95.0))/100.0
    all_isolates,_=load_matrix(run_dir,organism); selected=subset_samples(all_isolates,list(request.get("samples",[]))); retained={r["isolate_id"]:r for r in retained_rows(run_dir,organism)}; selected=[x for x in selected if x in retained]
    target_records,key=_gene_references(run_dir,organism,genes); target_ref=out/"target_references.fasta"; write_fasta(target_ref,target_records)
    flank_records={}; flank_meta={}
    for iso in selected:
        row=retained[iso]; assembly=Path(row.get("assembly","") or "missing"); locations=_locate(target_ref,assembly) if assembly.is_file() else {}
        if analyze_flanks and assembly.is_file():
            features=_gff_features(Path(row.get("gff","") or "missing"))
            for meta in key:
                rid=meta["reference_id"]
                for offset,feature in _neighbor_features(features,locations.get(rid,{}),flank_count) if rid in locations else []:
                    logical=f"{meta['Gene']}|{offset:+d}|{feature['name']}"; seq=_feature_sequence(assembly,feature)
                    if seq and logical not in flank_records:
                        flank_records[logical]=seq; flank_meta[logical]=(meta["Gene"],offset,str(feature["name"]),str(feature["locus_tag"]))
    records=list(target_records)
    for number,(logical,seq) in enumerate(flank_records.items(),len(records)+1):
        rid=f"CGF{number:08d}"; records.append((rid,seq)); parent,offset,name,locus=flank_meta[logical]; key.append({"reference_id":rid,"Gene":logical,"feature_type":"flank","parent_gene":parent,"flank_offset":str(offset),"annotation":name,"locus_tag":locus})
    reference=out/"variant_references.fasta"; write_fasta(reference,records); run(["bwa","index",str(reference)])
    refs=dict(records); output=[]; sequences={}
    for iso in selected:
        row=retained[iso]; sample_out=out/"evidence"/safe_name(iso); sample_out.mkdir(parents=True,exist_ok=True); bam=sample_out/"reads.bam"
        map_reads(reference,row["R1"],row["R2"],bam,threads,int(request.get("min_mapq",20)),sample_out/"bwa_mem.log")
        cov=coverage(bam,int(request.get("min_mapq",20))); fa=consensus(reference,bam,sample_out/"reads",float(request.get("min_depth",5)),int(request.get("min_mapq",20)),int(request.get("basequal",30))); cons=read_fasta(fa); aln=align_identity(reference,fa); var=_variant_counts(sample_out/"reads.vcf.gz")
        assembly=Path(row.get("assembly","") or "missing"); locations=_locate(reference,assembly) if assembly.is_file() else {}; features=_gff_features(Path(row.get("gff","") or "missing"));
        for meta in key:
            rid=meta["reference_id"]; c=cov.get(rid,{}); a=aln.get(rid)
            if not a and float(c.get("mapped_reads",0))>0: a=fixed_coordinate_identity(refs[rid],cons.get(rid,""))
            identity=None if not a else a.get("identity"); breadth=float(c.get("breadth",0)); location=locations.get(rid); neighbors=_neighbor_features(features,location,flank_count) if location else []
            vc=var.get(rid,{"snps":0,"mnps":0,"insertions":0,"deletions":0,"inserted_bases":0,"deleted_bases":0}); trunc=max(0,len(refs[rid])-int(c.get("covered_bases",0)))
            seq=cons.get(rid,""); digest=hashlib.sha256(seq.encode()).hexdigest() if seq and identity is not None and float(identity)>=min_similarity else ""
            if digest: sequences[(meta["Gene"],digest)]=seq
            output.append({"isolate_id":iso,"organism":organism,**meta,"variant_hash":digest,"percent_identity":"NA" if identity is None else float(identity)*100,"percent_coverage":breadth*100,"shared_bases":0 if not a else int(a.get("identical_positions",0)),"reference_length":len(refs[rid]),"snps":vc["snps"],"mnps":vc["mnps"],"insertions":vc["insertions"],"deletions":vc["deletions"],"inserted_bases":vc["inserted_bases"],"deleted_bases":vc["deleted_bases"],"truncation_bases":trunc,"location_status":"located" if location else "not_located" if assembly.is_file() else "no_assembly","contig":"" if not location else location["contig"],"start":"" if not location else location["start"],"end":"" if not location else location["end"],"strand":"" if not location else location["strand"],"flanking_genes":";".join(f"{offset:+d}:{feature['name']}" for offset,feature in neighbors)})
    variant_ids={}; unique_by_gene={}
    for gene,digest in sorted(sequences):
        index=unique_by_gene.get(gene,0)+1; unique_by_gene[gene]=index; variant_ids[(gene,digest)]=f"{safe_name(gene)}_V{index:03d}"
    for row in output: row["variant_id"]=variant_ids.get((row["Gene"],row.pop("variant_hash")),"")
    fields=["isolate_id","organism","Gene","feature_type","parent_gene","flank_offset","annotation","locus_tag","variant_id","percent_identity","percent_coverage","shared_bases","reference_length","snps","mnps","insertions","deletions","inserted_bases","deleted_bases","truncation_bases","location_status","contig","start","end","strand","flanking_genes"]
    write_tsv(out/"gene_variants.tsv",fields,output)
    counts=Counter(r["variant_id"] for r in output if r["variant_id"]); write_tsv(out/"variant_summary.tsv",["variant_id","Gene","n_isolates"],[[vid,next(r["Gene"] for r in output if r["variant_id"]==vid),n] for vid,n in sorted(counts.items())])
    plot_records=[(variant_ids[key],seq) for key,seq in sequences.items()]; _plot_variant_alignment(out,plot_records,counts)

CLASSIC=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f"]
MUTED=["#4878a8","#d9874c","#6a9f58","#c65f5f","#8b6f9f","#9c755f","#b279a2","#79706e"]

def _colors(labels: list[str], scheme: str, custom: dict[str,str]) -> dict[str,str]:
    palette=CLASSIC if scheme=="classic" else MUTED
    return {label:custom.get(label,palette[i%len(palette)]) for i,label in enumerate(labels)}

def _write_itol_binary(path: Path, label: str, fields: list[str], rows: list[tuple[str,list[int]]], colors: dict[str,str]) -> None:
    with path.open("w") as h:
        h.write("DATASET_BINARY\nSEPARATOR\tTAB\nDATASET_LABEL\t"+label+"\nCOLOR\t#333333\n")
        h.write("FIELD_LABELS\t"+"\t".join(fields)+"\nFIELD_COLORS\t"+"\t".join(colors[x] for x in fields)+"\nFIELD_SHAPES\t"+"\t".join("1" for _ in fields)+"\nDATA\n")
        for isolate,values in rows: h.write(isolate+"\t"+"\t".join(map(str,values))+"\n")

def _write_itol_strip(path: Path, label: str, values: dict[str,str], colors: dict[str,str]) -> None:
    with path.open("w") as h:
        h.write(f"DATASET_COLORSTRIP\nSEPARATOR\tTAB\nDATASET_LABEL\t{label}\nCOLOR\t#333333\nLEGEND_TITLE\t{label}\n")
        h.write("LEGEND_SHAPES\t"+"\t".join("1" for _ in colors)+"\nLEGEND_COLORS\t"+"\t".join(colors.values())+"\nLEGEND_LABELS\t"+"\t".join(colors)+"\nDATA\n")
        for isolate,value in values.items(): h.write(f"{isolate}\t{colors[value]}\t{value}\n")

def make_itol(request: dict[str,object]) -> None:
    run_dir=Path(str(request["run_dir"])); out=Path(str(request["output_dir"])); organism=str(request["organism"]); genes=list(request.get("genes",[])); scheme=str(request.get("color_scheme","classic")); custom=dict(request.get("custom_colors",{}))
    isolates,matrix=load_matrix(run_dir,organism); isolates=subset_samples(isolates,list(request.get("samples",[])))
    if genes:
        genes=select_genes(matrix,genes); colors=_colors(genes,scheme,custom); _write_itol_binary(out/"itol_gene_presence_absence.txt","Gene presence absence",genes,[(i,[matrix[g][i] for g in genes]) for i in isolates],colors)
    for key,column,group_column,label,filename in (("operon","operon_id","operon_name","Operon type","itol_operon_types.txt"),("variants","variant_id","Gene","Variant type","itol_variant_types.txt")):
        source=str(request.get(key,""))
        if not source: continue
        rows=read_tsv(Path(source)); groups=sorted({r.get(group_column,"") for r in rows})
        for group in groups:
            values={r["isolate_id"]:r[column] for r in rows if r.get(group_column,"")==group and r.get("isolate_id") in isolates and r.get(column)}
            if not values: continue
            labels=sorted(set(values.values())); target=filename if len(groups)==1 else f"{Path(filename).stem}.{safe_name(group)}.txt"
            _write_itol_strip(out/target,f"{label}: {group}",values,_colors(labels,scheme,custom))

def run_request(request_path: Path) -> None:
    request=load_json(request_path); out=Path(str(request["output_dir"])); out.mkdir(parents=True,exist_ok=True)
    kind=str(request["utility"])
    if kind=="get_samples": get_samples(request)
    elif kind=="get_differential_genes": differential_genes(request)
    elif kind=="get_operon": get_operon(request)
    elif kind=="get_variants": get_variants(request)
    elif kind=="diagnose_call":
        from .diagnostics import diagnose_call
        diagnose_call(request)
    elif kind=="itol": make_itol(request)
    else: raise SystemExit(f"Unknown CleanGene utility: {kind}")
