from __future__ import annotations
import subprocess, json
from pathlib import Path
from .fasta import read_fasta, write_fasta
from .util import read_tsv, run, write_tsv

EVIDENCE_VERSION="2"

METRIC_FIELDS=["evidence_version","family_breadth","reconstructed_coverage","reference_id","Gene","initial_call","validated_call","evidence_state","validation_state","decision_reason","decision_metrics","junction_identity","junction_spanning_alignment","flank_anchor_length","sequence_resolution","final_call_source","breadth","percent_coverage","mean_depth","normalized_depth","identity","percent_identity","identity_method","reconstructed_length","identical_positions","aligned_positions","reference_length","orf_integrity","mapped_reads","unique_mapped_reads","ambiguous_mapped_reads","mean_mapping_quality","assembly_scaffold","cds_start","cds_end","cds_strand","contig_edge","left_flank_locus","right_flank_locus","arbitration_status","arbitration_reason"]

def map_reads(reference: Path,r1: str,r2: str,bam: Path,threads: int,min_mapq: int,log: Path,*,retain_ambiguous: bool=False) -> None:
    bam.parent.mkdir(parents=True,exist_ok=True); log.parent.mkdir(parents=True,exist_ok=True)
    marker=bam.with_suffix(".mapping.json")
    signature={"version":EVIDENCE_VERSION,"min_mapq":min_mapq,"retain_ambiguous":retain_ambiguous,"inputs":[[str(Path(x).resolve()),Path(x).stat().st_size,Path(x).stat().st_mtime_ns] for x in (reference,r1,r2)]}
    if marker.is_file() and bam.is_file() and Path(str(bam)+".bai").is_file():
        try:
            if json.loads(marker.read_text())==signature: return
        except (ValueError,OSError): pass
    marker.unlink(missing_ok=True)
    view=["samtools","view","-u"]
    if not retain_ambiguous: view += ["-F","3332","-q",str(min_mapq)]
    view.append("-")
    with log.open("w") as err:
        bwa_command=["bwa","mem"]+(["-a"] if retain_ambiguous else [])+["-t",str(threads),str(reference),r1,r2]
        bwa=subprocess.Popen(bwa_command,stdout=subprocess.PIPE,stderr=err)
        sam=subprocess.Popen(view,stdin=bwa.stdout,stdout=subprocess.PIPE,stderr=err)
        assert bwa.stdout is not None and sam.stdout is not None; bwa.stdout.close()
        sort=subprocess.Popen(["samtools","sort","-@",str(max(1,threads-1)),"-o",str(bam),"-"],stdin=sam.stdout,stderr=err)
        sam.stdout.close(); statuses=(sort.wait(),sam.wait(),bwa.wait())
    if any(statuses): raise subprocess.CalledProcessError(next(x for x in statuses if x),"bwa mem | samtools view | samtools sort")
    run(["samtools","index",str(bam)])
    temporary=marker.with_suffix(".tmp"); temporary.write_text(json.dumps(signature)); temporary.replace(marker)

def coverage(bam: Path,min_mapq: int,*,include_ambiguous: bool=False) -> dict[str,dict[str,float]]:
    p=subprocess.run(["samtools","coverage","-q",str(0 if include_ambiguous else min_mapq),"--ff",str(1540 if include_ambiguous else 3844),str(bam)],check=True,capture_output=True,text=True); out={}
    for line in p.stdout.splitlines():
        if not line or line.startswith("#"): continue
        f=line.split("\t")
        if len(f)>=9: out[f[0]]={"mapped_reads":float(f[3]),"covered_bases":float(f[4]),"breadth":float(f[5])/100,"mean_depth":float(f[6]),"mean_base_quality":float(f[7]),"mean_mapping_quality":float(f[8])}
    return out

def region_coverage(bam: Path,contig: str,start: int,end: int,min_mapq: int) -> dict[str,float]:
    region=f"{contig}:{start}-{end}"; length=max(1,end-start+1)
    p=subprocess.run(["samtools","depth","-aa","-d","0","-G","2048","-Q",str(min_mapq),"-r",region,str(bam)],check=True,capture_output=True,text=True)
    depths=[int(line.rsplit("\t",1)[1]) for line in p.stdout.splitlines() if line]; depths += [0]*max(0,length-len(depths))
    mapped=int(subprocess.run(["samtools","view","-c","-F","3844","-q",str(min_mapq),str(bam),region],check=True,capture_output=True,text=True).stdout or 0)
    return {"mapped_reads":mapped,"breadth":sum(x>0 for x in depths)/length,"mean_depth":sum(depths)/length}

def representative_depth(cov: dict[str,dict[str,float]],assembly: Path) -> float:
    values=[(float(cov.get(k,{}).get("mean_depth",0)),len(v)) for k,v in read_fasta(assembly).items()]; total=sum(n for _,n in values); seen=0
    for depth,length in sorted(values):
        seen+=length
        if seen>=total/2: return depth
    return 0.0

def low_depth_bed(bam: Path,path: Path,min_depth: float,min_mapq: int) -> None:
    p=subprocess.run(["samtools","depth","-aa","-d","0","-G","2048","-Q",str(min_mapq),str(bam)],check=True,capture_output=True,text=True); intervals=[]; cur=None
    for line in p.stdout.splitlines():
        chrom,pos,dep=line.split("\t")[:3]; pos0=int(pos)-1
        if float(dep)>=min_depth:
            if cur: intervals.append(cur); cur=None
        elif cur and cur[0]==chrom and cur[2]==pos0: cur=(chrom,cur[1],pos0+1)
        else:
            if cur: intervals.append(cur)
            cur=(chrom,pos0,pos0+1)
    if cur: intervals.append(cur)
    with path.open("w") as handle:
        for x in intervals: handle.write(f"{x[0]}\t{x[1]}\t{x[2]}\n")

def consensus(reference: Path,bam: Path,prefix: Path,min_depth: float,min_mapq: int,basequal: int) -> Path:
    vcf=prefix.with_suffix(".vcf.gz"); mask=prefix.with_suffix(".low_depth.bed"); fa=prefix.with_suffix(".consensus.fasta")
    mp=subprocess.Popen(["bcftools","mpileup","-Ou","-f",str(reference),"-q",str(min_mapq),"-Q",str(basequal),"-d","100000",str(bam)],stdout=subprocess.PIPE); assert mp.stdout is not None
    call=subprocess.run(["bcftools","call","-mv","--ploidy","1","-Oz","-o",str(vcf)],stdin=mp.stdout); mp.stdout.close(); status=mp.wait()
    if status or call.returncode: raise subprocess.CalledProcessError(status or call.returncode,"bcftools")
    run(["bcftools","index","--force",str(vcf)]); low_depth_bed(bam,mask,min_depth,min_mapq)
    cmd=["bcftools","consensus","-c",str(prefix.with_suffix(".chain")),"-f",str(reference)]+(["-m",str(mask)] if mask.stat().st_size else [])+[str(vcf)]
    with fa.open("w") as out: subprocess.run(cmd,check=True,stdout=out)
    return fa

def fixed_coordinate_identity(ref: str,cons: str) -> dict[str,object]|None:
    ref=ref.upper(); cons=cons.upper(); compared=same=0
    for a,b in zip(ref,cons):
        if a=="N" or b=="N": continue
        compared+=1; same+=int(a==b)
    compared+=abs(len(ref)-len(cons))
    return None if not compared else {"identity":same/compared,"identical_positions":same,"aligned_positions":compared,"identity_method":"fixed_coordinate_global"}

def align_identity(reference: Path,consensus_fa: Path) -> dict[str,dict[str,float]]:
    """Compatibility helper used by downstream variant analyses."""
    p=subprocess.run(["minimap2","-x","asm5","--secondary=no","-c",str(reference),str(consensus_fa)],check=True,capture_output=True,text=True); best={}
    for line in p.stdout.splitlines():
        f=line.split("\t")
        if len(f)<12 or f[0]!=f[5]: continue
        matches,block=int(f[9]),int(f[10]); candidate={"identity":matches/block if block else None,"identical_positions":matches,"aligned_positions":block,"consensus_mapq":int(f[11]),"identity_method":"minimap2_asm5"}
        if f[5] not in best or matches>best[f[5]]["identical_positions"]: best[f[5]]=candidate
    return best

def orf_integrity(seq: str) -> str:
    seq=seq.upper()
    if not seq or any(base not in "ACGT" for base in seq): return "unresolved"
    if not seq or len(seq)%3: return "disrupted"
    return "intact" if all(seq[i:i+3] not in {"TAA","TAG","TGA"} for i in range(0,len(seq)-3,3)) else "internal_stop"

def best_sequence_match(reference_seq: str, contigs: Path, work: Path) -> dict[str,object]|None:
    work.mkdir(parents=True,exist_ok=True)
    reference=work/"target.fasta"; write_fasta(reference,[("target",reference_seq)])
    p=subprocess.run(["minimap2","-x","asm5","--secondary=no","-c",str(reference),str(contigs)],check=True,capture_output=True,text=True); best=None
    for line in p.stdout.splitlines():
        f=line.split("\t")
        if len(f)<12: continue
        candidate={"identity":int(f[9])/int(f[10]) if int(f[10]) else 0,"aligned_length":int(f[10]),"identical_positions":int(f[9]),"query_aligned_length":int(f[3])-int(f[2]),"reference_length":len(reference_seq),"breadth":min(1.0,(int(f[8])-int(f[7]))/max(1,len(reference_seq))),"contig":f[0],"reference_start":int(f[7]),"reference_end":int(f[8]),"alignment_cigar":next((x[5:] for x in f[12:] if x.startswith("cg:Z:")),"")}
        if best is None or (candidate["breadth"],candidate["identity"])>(best["breadth"],best["identity"]): best=candidate
    return best

def targeted_local_reconstruction(*,bam: Path,region: str,reference_seq: str,outdir: Path,threads: int,flank_junction: str="",junction_offset: int=0,max_reads: int=100000,memory_gb: int=12,deletion_identity: float=.95,deletion_anchor: int=50) -> dict[str,object]:
    """Recruit alignments and their mates, locally assemble, and resolve a target."""
    outdir.mkdir(parents=True,exist_ok=True); names=outdir/"read_names.txt"
    command=["samtools","view",str(bam),region]; read_names=set()
    with subprocess.Popen(command,stdout=subprocess.PIPE,text=True) as process:
        assert process.stdout is not None
        for line in process.stdout:
            if line: read_names.add(line.split("\t",1)[0])
            if len(read_names)>max_reads:
                process.terminate(); process.wait()
                return {"status":"recruitment_limit"}
        status=process.wait()
    if status: raise subprocess.CalledProcessError(status,command)
    if not read_names: return {"status":"no_recruited_reads"}
    read_names=sorted(read_names)
    names.write_text("".join(f"{name}\n" for name in read_names))
    recruited=outdir/"recruited.bam"; run(["samtools","view","-b","-F","2304","-N",str(names),"-o",str(recruited),str(bam)])
    collated=outdir/"recruited.collated.bam"
    run(["samtools","collate","-o",str(collated),str(recruited)])
    r1=outdir/"recruited_R1.fastq"; r2=outdir/"recruited_R2.fastq"; singles=outdir/"recruited_singletons.fastq"; other=outdir/"recruited_other.fastq"
    run(["samtools","fastq","-n","-1",str(r1),"-2",str(r2),"-0",str(other),"-s",str(singles),str(collated)],stdout=outdir/"samtools_fastq.stdout",stderr=outdir/"samtools_fastq.stderr")
    assembly=outdir/"spades"; command=["spades.py","--only-assembler","--careful","-m",str(memory_gb),"-t",str(threads),"-o",str(assembly)]
    has_reads=False
    if r1.stat().st_size and r2.stat().st_size: command += ["-1",str(r1),"-2",str(r2)]; has_reads=True
    if singles.stat().st_size or other.stat().st_size:
        single_input=outdir/"all_singletons.fastq"
        with single_input.open("wb") as handle:
            for source in (singles,other):
                with source.open("rb") as source_handle:
                    import shutil
                    shutil.copyfileobj(source_handle,handle)
        command += ["-s",str(single_input)]; has_reads=True
    if not has_reads: return {"status":"no_reconstructed_reads"}
    run(command,stdout=outdir/"spades.stdout",stderr=outdir/"spades.stderr"); contigs=assembly/"contigs.fasta"
    if not contigs.is_file() or not read_fasta(contigs): return {"status":"no_contigs"}
    match=best_sequence_match(reference_seq,contigs,outdir/"candidate_match") if reference_seq else None
    deletion=best_sequence_match(flank_junction,contigs,outdir/"deletion_match") if flank_junction else None
    if deletion: deletion["flank_anchor_length"]=junction_anchor_length(deletion,junction_offset)
    deletion_spanned=bool(deletion and junction_offset>=deletion_anchor and len(flank_junction)-junction_offset>=deletion_anchor and deletion["identity"]>=deletion_identity and spans_junction(deletion,junction_offset,deletion_anchor))
    return {"status":"reconstructed","candidate":match,"deletion_spanned":deletion_spanned,"deletion":deletion}

def classify_gene_evidence(*,initial_call: int=0,mapped_reads: float,breadth: float,mean_depth: float,identity: float|None,min_breadth: float=.95,min_depth: float=5,min_identity: float=.95,truncation_breadth: float=.70,divergent_breadth: float=.90,divergent_identity: float=.90,unique_reads: float|None=None,ambiguous_reads: float=0) -> dict[str,object]:
    if ambiguous_reads>0 and (unique_reads or 0)==0: state,call,source="ambiguous_multimap","","arbitration_pending"
    elif not mapped_reads or not breadth: state,call,source="not_detected","" if initial_call else 0,"arbitration_pending" if initial_call else "read_validation"
    elif mean_depth<min_depth or identity is None: state,call,source="insufficient_evidence","","initial_call_unresolved"
    elif breadth>=min_breadth and identity>=min_identity: state,call,source="confirmed_present",1,"own_locus_read_validation" if initial_call else "pangenome_read_recovery"
    elif breadth>=divergent_breadth and divergent_identity<=identity<min_identity: state,call,source="divergent_variant",1,"read_validation"
    elif breadth>=truncation_breadth and identity>=min_identity: state,call,source="possible_truncation",1 if initial_call else "","arbitration_pending"
    else: state,call,source="partial_homolog",0,"arbitration_pending" if initial_call else "read_validation"
    from .validation_summary import STATE_METRICS
    metrics=STATE_METRICS[state]
    if state=="not_detected": metrics=("mapped_reads",) if not mapped_reads else ("breadth",)
    elif state=="insufficient_evidence": metrics=("mean_depth",) if mean_depth<min_depth else ("identity",)
    return {"evidence_state":state,"validation_state":state,"validated_call":call,"final_call_source":source,"decision_reason":state.replace("_"," "),"decision_metrics":";".join(metrics)}

def validation_decision_logic_rows(min_breadth="0.95",min_depth="5",min_identity="0.95") -> list[list[str]]:
    return [["confirmed_present",f"breadth >= {min_breadth}; identity >= {min_identity}; depth >= {min_depth}","1","Intact sequence supported"],["possible_truncation","breadth 0.70 to confirmed threshold; high identity","initial positive: 1 pending arbitration","Possible endpoint or assembly break"],["divergent_variant","breadth >= 0.90; identity 0.90 to confirmed threshold","1","Divergent full-length allele"],["partial_homolog","breadth below 0.70 or weak similarity","0","Related sequence, not an intact gene"],["ambiguous_multimap","family mappings but no unique assignment","preserve/arbitrate","Family present; exact cluster unresolved"],["not_detected","no meaningful read evidence","initial positive: arbitrate; otherwise 0","Absence not proven without locus evidence"],["confirmed_absent_locus","flank reconstruction spans deletion","0","Physical deletion junction supported"]]

def _slice(seqs: dict[str,str],row: dict[str,str]) -> str:
    seq=seqs.get(row.get("assembly_scaffold",""),"")[int(row.get("cds_start") or 1)-1:int(row.get("cds_end") or 0)]
    return seq.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1] if row.get("cds_strand")=="-" else seq

def validate_isolate(reference: Path,key_tsv: Path,locus_tsv: Path,assembly: Path,r1: str,r2: str,outdir: Path,threads: int,min_breadth: float,min_depth: float,min_identity: float,min_mapq: int,basequal: int,*,initial_calls: dict[str,int]|None=None,truncation_breadth: float=.70,divergent_breadth: float=.90,divergent_identity: float=.90) -> None:
    outdir.mkdir(parents=True,exist_ok=True); keys=read_tsv(key_tsv); locus_rows=read_tsv(locus_tsv); loci={}; initial_calls=initial_calls or {}
    for locus_row in locus_rows: loci.setdefault(locus_row["Gene"],[]).append(locus_row)
    if assembly.suffix==".gz":
        uncompressed=outdir/"own_assembly.fasta"
        if not uncompressed.is_file(): write_fasta(uncompressed,list(read_fasta(assembly).items()))
        assembly=uncompressed
    if not Path(str(assembly)+".bwt").is_file(): run(["bwa","index",str(assembly)],stdout=outdir/"own_assembly_bwa_index.stdout",stderr=outdir/"own_assembly_bwa_index.stderr")
    if not Path(str(assembly)+".fai").is_file(): run(["samtools","faidx",str(assembly)])
    own=outdir/"own_assembly_reads.bam"; map_reads(assembly,r1,r2,own,threads,min_mapq,outdir/"own_assembly_bwa.log",retain_ambiguous=True); own_cov=coverage(own,min_mapq); chrom_depth=representative_depth(own_cov,assembly)
    own_cons=read_fasta(consensus(assembly,own,outdir/"own_assembly_reads",min_depth,min_mapq,basequal)); assembly_seqs=read_fasta(assembly)
    locus_metrics=locus_coverage(own,locus_rows,outdir,min_mapq)
    def locus_key(row): return (row["assembly_scaffold"],int(row["cds_start"]),int(row["cds_end"]))
    search=outdir/"pangenome_reads.bam"; map_reads(reference,r1,r2,search,threads,min_mapq,outdir/"pangenome_bwa.log",retain_ambiguous=True); search_cov=coverage(search,min_mapq,include_ambiguous=True); search_unique_cov=coverage(search,min_mapq); search_cons=read_fasta(consensus(reference,search,outdir/"pangenome_reads",min_depth,min_mapq,basequal)); refs=read_fasta(reference); rows=[]
    for key in keys:
        gene=key["Gene"]; initial=int(initial_calls.get(gene,key.get("initial_call",0))); candidates=loci.get(gene,[]) if initial else []; locus=max(candidates,key=lambda r:(locus_metrics[locus_key(r)]["breadth"],locus_metrics[locus_key(r)]["mean_depth"])) if candidates else None
        if locus:
            c=locus_metrics[locus_key(locus)]; refseq=_slice(assembly_seqs,locus); reconstructed=consensus_locus(own_cons,locus,outdir/"own_assembly_reads.chain"); ident=sequence_identity(refseq,reconstructed); unique=int(c["mapped_reads"]); ambiguous=int(c["ambiguous_mapped_reads"]); resolution="exact"
        else:
            c=search_cov.get(key["reference_id"],{}); refseq=refs.get(key["reference_id"],""); reconstructed=search_cons.get(key["reference_id"],"") if c.get("breadth",0) else ""; ident=sequence_identity(refseq,reconstructed) if reconstructed else None; total=int(c.get("mapped_reads",0)); unique=int(search_unique_cov.get(key["reference_id"],{}).get("mapped_reads",0)); ambiguous=max(0,total-unique); resolution="reconstructed" if unique else "family_only" if ambiguous else "unresolved"
        family_breadth=float(c.get("breadth",0))
        if not locus and unique:
            c=search_unique_cov.get(key["reference_id"],{})
        identity=None if not ident else ident["identity"]; decision=classify_gene_evidence(initial_call=initial,mapped_reads=float(c.get("mapped_reads",0)),breadth=float(c.get("breadth",0)),mean_depth=float(c.get("mean_depth",0)),identity=identity,min_breadth=min_breadth,min_depth=min_depth,min_identity=min_identity,truncation_breadth=truncation_breadth,divergent_breadth=divergent_breadth,divergent_identity=divergent_identity,unique_reads=unique,ambiguous_reads=ambiguous)
        if ambiguous and float(c.get("breadth",0))<min_breadth and family_breadth>=min_breadth:
            decision.update(evidence_state="ambiguous_multimap",validation_state="ambiguous_multimap",validated_call="",final_call_source="arbitration_pending",decision_reason="ambiguous multimap",decision_metrics="ambiguous_mapped_reads;breadth;family_breadth"); resolution="family_only"
        if decision["evidence_state"]=="ambiguous_multimap": resolution="family_only"
        elif decision["evidence_state"] in {"not_detected","insufficient_evidence"}: resolution="unresolved"
        row={f:"" for f in METRIC_FIELDS}; row.update(key); row.update(decision); row.update({"evidence_version":EVIDENCE_VERSION,"family_breadth":family_breadth,"reconstructed_coverage":min(1.,sum(b in "ACGTacgt" for b in reconstructed)/len(refseq)) if refseq else 0.,"initial_call":initial,"sequence_resolution":resolution,"breadth":c.get("breadth",0),"percent_coverage":float(c.get("breadth",0))*100,"mean_depth":c.get("mean_depth",0),"normalized_depth":float(c.get("mean_depth",0))/chrom_depth if chrom_depth else "","identity":"NA" if identity is None else identity,"percent_identity":"NA" if identity is None else identity*100,"identity_method":ident.get("identity_method","") if ident else "","reconstructed_length":len(reconstructed.replace("N","")),"identical_positions":ident.get("identical_positions","") if ident else "","aligned_positions":ident.get("aligned_positions","") if ident else "","reference_length":len(refseq),"orf_integrity":orf_integrity(reconstructed),"mapped_reads":unique+ambiguous,"unique_mapped_reads":unique,"ambiguous_mapped_reads":ambiguous,"mean_mapping_quality":c.get("mean_mapping_quality","")});
        if locus: row.update({k:locus.get(k,"") for k in ("assembly_scaffold","cds_start","cds_end","cds_strand","contig_edge","left_flank_locus","right_flank_locus")})
        row["arbitration_status"]="pending" if row["final_call_source"]=="arbitration_pending" or (not initial and row["validated_call"]==1) or (initial and row["evidence_state"]=="insufficient_evidence") else "not_required"; rows.append(row)
    write_tsv(outdir/"metrics.tsv",METRIC_FIELDS,rows)


def sequence_identity(ref: str, reconstructed: str) -> dict[str,object]|None:
    """Compare reconstructed bases with gaps, excluding unknown reference positions."""
    ref=ref.upper(); reconstructed=reconstructed.upper()
    if not reconstructed or not any(b in "ACGT" for b in reconstructed): return None
    # Strip identical ends before dynamic programming; most own-locus calls
    # therefore require no alignment matrix. Bound exceptional divergent cases.
    left=0; right=0; limit=min(len(ref),len(reconstructed))
    while left<limit and ref[left]==reconstructed[left]: left+=1
    while right<limit-left and ref[len(ref)-right-1]==reconstructed[len(reconstructed)-right-1]: right+=1
    ends=ref[:left]+(ref[len(ref)-right:] if right else "")
    a=ref[left:len(ref)-right if right else len(ref)]
    b=reconstructed[left:len(reconstructed)-right if right else len(reconstructed)]
    if len(a)*len(b)>2_000_000: return None
    # Entries contain edit cost, matches, and observed alignment columns.
    previous=[(j,0,sum(x in "ACGT" for x in b[:j])) for j in range(len(b)+1)]
    for i,x in enumerate(a,1):
        current=[(i,0,sum(t in "ACGT" for t in a[:i]))]
        for j,y in enumerate(b,1):
            known=x in "ACGT" and y in "ACGT"; equal=x==y and known
            cost,matches,columns=previous[j-1]
            diagonal=(cost+int(known and not equal),matches+int(equal),columns+int(known))
            cost,matches,columns=previous[j]; deletion=(cost+1,matches,columns+int(x in "ACGT"))
            cost,matches,columns=current[j-1]; insertion=(cost+1,matches,columns+int(y in "ACGT"))
            current.append(min((diagonal,deletion,insertion),key=lambda entry:(entry[0],-entry[1],entry[2])))
        previous=current
    _,same,compared=previous[-1]; end_count=sum(x in "ACGT" for x in ends); same+=end_count; compared+=end_count
    return None if not compared else {"identity":same/compared,"identical_positions":same,"aligned_positions":compared,"identity_method":"reconstructed_global_edit_alignment"}


def consensus_locus(seqs: dict[str,str], locus: dict[str,str], chain: Path) -> str:
    """Lift reference CDS endpoints through the bcftools consensus chain."""
    start=int(locus["cds_start"])-1; end=int(locus["cds_end"]); chrom=locus["assembly_scaffold"]
    blocks=[]; active=False; target=query=0
    for line in chain.read_text().splitlines():
        fields=line.split()
        if not fields: continue
        if fields[0]=="chain":
            active=fields[2]==chrom; target=int(fields[5]); query=int(fields[10]); continue
        if not active: continue
        size=int(fields[0]); blocks.append((target,target+size,query,query+size))
        target+=size; query+=size
        if len(fields)==3:
            dt,dq=map(int,fields[1:]); blocks.append((target,target+dt,query,query+dq)); target+=dt; query+=dq
    def lift(pos: int) -> int:
        for a,b,c,d in blocks:
            if a<=pos<=b: return c+min(pos-a,d-c)
        raise ValueError(f"CDS coordinate {chrom}:{pos} is outside consensus chain")
    lifted={**locus,"cds_start":str(lift(start)+1),"cds_end":str(lift(end))}
    return _slice(seqs,lifted)


def junction_anchor_length(match: dict[str,object], offset: int) -> int:
    """Observed minimum aligned flank length in a contiguous junction block."""
    import re
    pos=int(match["reference_start"]); anchor=0
    for length,op in re.findall(r"(\d+)([MIDNSHP=X])",str(match.get("alignment_cigar",""))):
        length=int(length)
        if op in "M=X" and pos<=offset<=pos+length:
            anchor=max(anchor,min(offset-pos,pos+length-offset))
        if op in "MDN=X": pos+=length
    return anchor


def spans_junction(match: dict[str,object], offset: int, anchor: int=50) -> bool:
    """Require a contiguous aligned block on both sides; gaps cannot prove a join."""
    observed=junction_anchor_length(match,offset)
    return observed>0 and observed>=anchor


def locus_coverage(bam: Path,loci: list[dict[str,str]],outdir: Path,min_mapq: int) -> dict[tuple,dict[str,float]]:
    """Measure all CDSs in two streamed BAM passes, including overlapping CDSs."""
    import bisect, re
    regions={}
    for row in loci:
        key=(row['assembly_scaffold'],int(row['cds_start']),int(row['cds_end']))
        regions[key]={"mapped_reads":0,"ambiguous_mapped_reads":0,"covered":0,"depth_sum":0}
    if not regions: return {}
    by_contig={}
    for key in sorted(regions): by_contig.setdefault(key[0],[]).append(key)
    indices={}
    for chrom,keys in by_contig.items():
        maxima=[]; maximum=0
        for _,start,end in keys: maximum=max(maximum,end); maxima.append(maximum)
        indices[chrom]=([k[1] for k in keys],maxima,keys)
    def overlapping(chrom,start,end):
        if chrom not in indices: return
        starts,maxima,keys=indices[chrom]; i=bisect.bisect_right(starts,end)-1
        while i>=0 and maxima[i]>=start:
            if keys[i][2]>=start: yield keys[i]
            i-=1
    bed=outdir/'cds_regions.bed'
    bed.write_text(''.join(f'{chrom}\t{start-1}\t{end}\n' for chrom,start,end in regions))
    commands=[['samtools','depth','-d','0','-G','2048','-Q',str(min_mapq),'-b',str(bed),str(bam)],['samtools','view','-M','-L',str(bed),'-F','1540',str(bam)]]
    for mode,command in enumerate(commands):
        with subprocess.Popen(command,stdout=subprocess.PIPE,text=True) as process:
            assert process.stdout is not None
            for line in process.stdout:
                fields=line.rstrip().split('\t')
                if mode==0:
                    chrom,pos,depth=fields[:3]; pos=int(pos); depth=int(depth)
                    for key in overlapping(chrom,pos,pos): regions[key]['covered']+=int(depth>0); regions[key]['depth_sum']+=depth
                else:
                    flag=int(fields[1]); chrom=fields[2]; start=int(fields[3]); mapq=int(fields[4])
                    span=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',fields[5]) if op in 'MDN=X')
                    kind='mapped_reads' if mapq>=min_mapq and not flag&2304 else 'ambiguous_mapped_reads'
                    for key in overlapping(chrom,start,start+span-1): regions[key][kind]+=1
            status=process.wait()
        if status: raise subprocess.CalledProcessError(status,command)
    for (_,start,end),values in regions.items():
        values['breadth']=values.pop('covered')/(end-start+1); values['mean_depth']=values.pop('depth_sum')/(end-start+1)
    return regions
