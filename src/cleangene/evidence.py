from __future__ import annotations
import csv, subprocess
from pathlib import Path
from .fasta import read_fasta, write_fasta
from .util import run, write_tsv

METRIC_FIELDS=["reference_id","Gene","initial_call","validated_call","evidence_state","validation_state","decision_reason","sequence_resolution","final_call_source","breadth","percent_coverage","mean_depth","normalized_depth","identity","percent_identity","identity_method","reconstructed_length","identical_positions","aligned_positions","reference_length","orf_integrity","mapped_reads","unique_mapped_reads","ambiguous_mapped_reads","mean_mapping_quality","assembly_scaffold","cds_start","cds_end","cds_strand","contig_edge","left_flank_locus","right_flank_locus","arbitration_status","arbitration_reason"]

def map_reads(reference: Path,r1: str,r2: str,bam: Path,threads: int,min_mapq: int,log: Path,*,retain_ambiguous: bool=False) -> None:
    bam.parent.mkdir(parents=True,exist_ok=True); log.parent.mkdir(parents=True,exist_ok=True)
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

def coverage(bam: Path,min_mapq: int) -> dict[str,dict[str,float]]:
    p=subprocess.run(["samtools","coverage",str(bam)],check=True,capture_output=True,text=True); out={}
    for line in p.stdout.splitlines():
        if not line or line.startswith("#"): continue
        f=line.split("\t")
        if len(f)>=9: out[f[0]]={"mapped_reads":float(f[3]),"covered_bases":float(f[4]),"breadth":float(f[5])/100,"mean_depth":float(f[6]),"mean_base_quality":float(f[7]),"mean_mapping_quality":float(f[8])}
    return out

def region_coverage(bam: Path,contig: str,start: int,end: int,min_mapq: int) -> dict[str,float]:
    region=f"{contig}:{start}-{end}"; length=max(1,end-start+1)
    p=subprocess.run(["samtools","depth","-aa","-d","0","-Q",str(min_mapq),"-r",region,str(bam)],check=True,capture_output=True,text=True)
    depths=[int(line.rsplit("\t",1)[1]) for line in p.stdout.splitlines() if line]; depths += [0]*max(0,length-len(depths))
    mapped=int(subprocess.run(["samtools","view","-c","-q",str(min_mapq),str(bam),region],check=True,capture_output=True,text=True).stdout or 0)
    return {"mapped_reads":mapped,"breadth":sum(x>0 for x in depths)/length,"mean_depth":sum(depths)/length}

def representative_depth(cov: dict[str,dict[str,float]],assembly: Path) -> float:
    values=[(float(cov.get(k,{}).get("mean_depth",0)),len(v)) for k,v in read_fasta(assembly).items()]; total=sum(n for _,n in values); seen=0
    for depth,length in sorted(values):
        seen+=length
        if seen>=total/2: return depth
    return 0.0

def low_depth_bed(bam: Path,path: Path,min_depth: float,min_mapq: int) -> None:
    p=subprocess.run(["samtools","depth","-aa","-d","0","-Q",str(min_mapq),str(bam)],check=True,capture_output=True,text=True); intervals=[]; cur=None
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
    cmd=["bcftools","consensus","-f",str(reference)]+(["-m",str(mask)] if mask.stat().st_size else [])+[str(vcf)]
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
    seq=seq.upper().replace("N","")
    if not seq or len(seq)%3: return "disrupted"
    return "intact" if all(seq[i:i+3] not in {"TAA","TAG","TGA"} for i in range(0,len(seq)-3,3)) else "internal_stop"

def best_sequence_match(reference_seq: str, contigs: Path, work: Path) -> dict[str,object]|None:
    reference=work/"target.fasta"; write_fasta(reference,[("target",reference_seq)])
    p=subprocess.run(["minimap2","-x","asm5","--secondary=no","-c",str(reference),str(contigs)],check=True,capture_output=True,text=True); best=None
    for line in p.stdout.splitlines():
        f=line.split("\t")
        if len(f)<12: continue
        candidate={"identity":int(f[9])/int(f[10]) if int(f[10]) else 0,"aligned_length":int(f[10]),"reference_length":len(reference_seq),"breadth":min(1.0,(int(f[8])-int(f[7]))/max(1,len(reference_seq))),"contig":f[0]}
        if best is None or (candidate["breadth"],candidate["identity"])>(best["breadth"],best["identity"]): best=candidate
    return best

def targeted_local_reconstruction(*,bam: Path,region: str,reference_seq: str,outdir: Path,threads: int,flank_junction: str="") -> dict[str,object]:
    """Recruit alignments and their mates, locally assemble, and resolve a target."""
    outdir.mkdir(parents=True,exist_ok=True); names=outdir/"read_names.txt"
    p=subprocess.run(["samtools","view",str(bam),region],check=True,capture_output=True,text=True)
    read_names=sorted({line.split("\t",1)[0] for line in p.stdout.splitlines() if line})
    if not read_names: return {"status":"no_recruited_reads"}
    names.write_text("".join(f"{name}\n" for name in read_names))
    recruited=outdir/"recruited.bam"; run(["samtools","view","-b","-F","2304","-N",str(names),"-o",str(recruited),str(bam)])
    r1=outdir/"recruited_R1.fastq"; r2=outdir/"recruited_R2.fastq"; singles=outdir/"recruited_singletons.fastq"; other=outdir/"recruited_other.fastq"
    run(["samtools","fastq","-n","-1",str(r1),"-2",str(r2),"-0",str(other),"-s",str(singles),str(recruited)],stdout=outdir/"samtools_fastq.stdout",stderr=outdir/"samtools_fastq.stderr")
    assembly=outdir/"spades"; command=["spades.py","--only-assembler","--careful","-t",str(threads),"-o",str(assembly)]
    has_reads=False
    if r1.stat().st_size and r2.stat().st_size: command += ["-1",str(r1),"-2",str(r2)]; has_reads=True
    if singles.stat().st_size or other.stat().st_size:
        single_input=singles if singles.stat().st_size else other; command += ["-s",str(single_input)]; has_reads=True
    if not has_reads: return {"status":"no_reconstructed_reads"}
    run(command,stdout=outdir/"spades.stdout",stderr=outdir/"spades.stderr"); contigs=assembly/"contigs.fasta"
    if not contigs.is_file() or not read_fasta(contigs): return {"status":"no_contigs"}
    match=best_sequence_match(reference_seq,contigs,outdir/"candidate_match") if reference_seq else None
    deletion=best_sequence_match(flank_junction,contigs,outdir/"deletion_match") if flank_junction else None
    deletion_spanned=bool(deletion and deletion["breadth"]>=.90 and deletion["identity"]>=.95)
    return {"status":"reconstructed","candidate":match,"deletion_spanned":deletion_spanned,"deletion":deletion}

def classify_gene_evidence(*,initial_call: int=0,mapped_reads: float,breadth: float,mean_depth: float,identity: float|None,min_breadth: float=.95,min_depth: float=5,min_identity: float=.95,truncation_breadth: float=.70,divergent_breadth: float=.90,divergent_identity: float=.90,unique_reads: float|None=None,ambiguous_reads: float=0) -> dict[str,object]:
    if ambiguous_reads>0 and (unique_reads or 0)==0: state,call,source="ambiguous_multimap","","arbitration_pending"
    elif not mapped_reads or not breadth: state,call,source="not_detected","" if initial_call else 0,"arbitration_pending" if initial_call else "read_validation"
    elif mean_depth<min_depth or identity is None: state,call,source="insufficient_evidence","","initial_call_unresolved"
    elif breadth>=min_breadth and identity>=min_identity: state,call,source="confirmed_present",1,"own_locus_read_validation" if initial_call else "pangenome_read_recovery"
    elif breadth>=divergent_breadth and identity>=divergent_identity: state,call,source="divergent_variant",1,"read_validation"
    elif breadth>=truncation_breadth and identity>=min_identity: state,call,source="possible_truncation",1 if initial_call else "","arbitration_pending"
    else: state,call,source="partial_homolog",0,"arbitration_pending" if initial_call else "read_validation"
    return {"evidence_state":state,"validation_state":state,"validated_call":call,"final_call_source":source,"decision_reason":state.replace("_"," ")}

def validation_decision_logic_rows(min_breadth="0.95",min_depth="5",min_identity="0.95") -> list[list[str]]:
    return [["confirmed_present",f"breadth >= {min_breadth}; identity >= {min_identity}; depth >= {min_depth}","1","Intact sequence supported"],["possible_truncation","breadth 0.70 to confirmed threshold; high identity","initial positive: 1 pending arbitration","Possible endpoint or assembly break"],["divergent_variant","breadth >= 0.90; identity 0.90 to confirmed threshold","1","Divergent full-length allele"],["partial_homolog","breadth below 0.70 or weak similarity","0","Related sequence, not an intact gene"],["ambiguous_multimap","family mappings but no unique assignment","preserve/arbitrate","Family present; exact cluster unresolved"],["not_detected","no meaningful read evidence","initial positive: arbitrate; otherwise 0","Absence not proven without locus evidence"],["confirmed_absent_locus","flank reconstruction spans deletion","0","Physical deletion junction supported"]]

def _slice(seqs: dict[str,str],row: dict[str,str]) -> str:
    seq=seqs.get(row.get("assembly_scaffold",""),"")[int(row.get("cds_start") or 1)-1:int(row.get("cds_end") or 0)]
    return seq.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1] if row.get("cds_strand")=="-" else seq

def validate_isolate(reference: Path,key_tsv: Path,locus_tsv: Path,assembly: Path,r1: str,r2: str,outdir: Path,threads: int,min_breadth: float,min_depth: float,min_identity: float,min_mapq: int,basequal: int,*,initial_calls: dict[str,int]|None=None,truncation_breadth: float=.70,divergent_breadth: float=.90,divergent_identity: float=.90) -> None:
    outdir.mkdir(parents=True,exist_ok=True); keys=list(csv.DictReader(key_tsv.open(newline=""),delimiter="\t")); loci={r["Gene"]:r for r in csv.DictReader(locus_tsv.open(newline=""),delimiter="\t")}; initial_calls=initial_calls or {}
    if assembly.suffix==".gz":
        uncompressed=outdir/"own_assembly.fasta"
        if not uncompressed.is_file(): write_fasta(uncompressed,list(read_fasta(assembly).items()))
        assembly=uncompressed
    if not Path(str(assembly)+".bwt").is_file(): run(["bwa","index",str(assembly)],stdout=outdir/"own_assembly_bwa_index.stdout",stderr=outdir/"own_assembly_bwa_index.stderr")
    if not Path(str(assembly)+".fai").is_file(): run(["samtools","faidx",str(assembly)])
    own=outdir/"own_assembly_reads.bam"; map_reads(assembly,r1,r2,own,threads,min_mapq,outdir/"own_assembly_bwa.log"); own_cov=coverage(own,min_mapq); chrom_depth=representative_depth(own_cov,assembly)
    own_cons=read_fasta(consensus(assembly,own,outdir/"own_assembly_reads",min_depth,min_mapq,basequal)); assembly_seqs=read_fasta(assembly)
    search=outdir/"pangenome_reads.bam"; map_reads(reference,r1,r2,search,threads,min_mapq,outdir/"pangenome_bwa.log",retain_ambiguous=True); search_cov=coverage(search,min_mapq); search_cons=read_fasta(consensus(reference,search,outdir/"pangenome_reads",min_depth,min_mapq,basequal)); refs=read_fasta(reference); rows=[]
    for key in keys:
        gene=key["Gene"]; initial=int(initial_calls.get(gene,key.get("initial_call",0))); locus=loci.get(gene) if initial else None
        if locus:
            c=region_coverage(own,locus["assembly_scaffold"],int(locus["cds_start"]),int(locus["cds_end"]),min_mapq); refseq=_slice(assembly_seqs,locus); reconstructed=_slice(own_cons,locus); ident=fixed_coordinate_identity(refseq,reconstructed); unique=int(c["mapped_reads"]); ambiguous=0; resolution="exact"
        else:
            c=search_cov.get(key["reference_id"],{}); refseq=refs.get(key["reference_id"],""); reconstructed=search_cons.get(key["reference_id"],"") if c.get("breadth",0) else ""; ident=fixed_coordinate_identity(refseq,reconstructed) if reconstructed else None; total=int(c.get("mapped_reads",0)); unique=int(subprocess.run(["samtools","view","-c","-F","2304","-q",str(min_mapq),str(search),key["reference_id"]],check=True,capture_output=True,text=True).stdout or 0); ambiguous=max(0,total-unique); resolution="reconstructed" if unique else "family_only" if ambiguous else "unresolved"
        identity=None if not ident else ident["identity"]; decision=classify_gene_evidence(initial_call=initial,mapped_reads=float(c.get("mapped_reads",0)),breadth=float(c.get("breadth",0)),mean_depth=float(c.get("mean_depth",0)),identity=identity,min_breadth=min_breadth,min_depth=min_depth,min_identity=min_identity,truncation_breadth=truncation_breadth,divergent_breadth=divergent_breadth,divergent_identity=divergent_identity,unique_reads=unique,ambiguous_reads=ambiguous)
        row={f:"" for f in METRIC_FIELDS}; row.update(key); row.update(decision); row.update({"initial_call":initial,"sequence_resolution":resolution,"breadth":c.get("breadth",0),"percent_coverage":float(c.get("breadth",0))*100,"mean_depth":c.get("mean_depth",0),"normalized_depth":float(c.get("mean_depth",0))/chrom_depth if chrom_depth else "","identity":"NA" if identity is None else identity,"percent_identity":"NA" if identity is None else identity*100,"identity_method":ident.get("identity_method","") if ident else "","reconstructed_length":len(reconstructed.replace("N","")),"identical_positions":ident.get("identical_positions","") if ident else "","aligned_positions":ident.get("aligned_positions","") if ident else "","reference_length":len(refseq),"orf_integrity":orf_integrity(reconstructed),"mapped_reads":int(c.get("mapped_reads",0)),"unique_mapped_reads":unique,"ambiguous_mapped_reads":ambiguous,"mean_mapping_quality":c.get("mean_mapping_quality","")});
        if locus: row.update({k:locus.get(k,"") for k in ("assembly_scaffold","cds_start","cds_end","cds_strand","contig_edge","left_flank_locus","right_flank_locus")})
        row["arbitration_status"]="pending" if row["final_call_source"]=="arbitration_pending" else "not_required"; rows.append(row)
    write_tsv(outdir/"metrics.tsv",METRIC_FIELDS,rows)
