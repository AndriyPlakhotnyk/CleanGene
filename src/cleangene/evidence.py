from __future__ import annotations
import csv, subprocess
from pathlib import Path
from .fasta import read_fasta
from .util import run, write_tsv

def map_reads(reference: Path, r1: str, r2: str, bam: Path, threads: int, min_mapq: int, log: Path) -> None:
    bam.parent.mkdir(parents=True, exist_ok=True); log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as err:
        bwa = subprocess.Popen(["bwa","mem","-t",str(threads),str(reference),r1,r2], stdout=subprocess.PIPE, stderr=err)
        view = subprocess.Popen(["samtools","view","-u","-F","3332","-q",str(min_mapq),"-"], stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=err)
        assert bwa.stdout is not None and view.stdout is not None; bwa.stdout.close()
        sort = subprocess.Popen(["samtools","sort","-@",str(max(1,threads-1)),"-o",str(bam),"-"], stdin=view.stdout, stderr=err)
        view.stdout.close(); s, v, b = sort.wait(), view.wait(), bwa.wait()
    if s or v or b: raise subprocess.CalledProcessError(s or v or b, "bwa mem | samtools view | samtools sort")
    run(["samtools","index",str(bam)])

def coverage(bam: Path, min_mapq: int) -> dict[str, dict[str, float]]:
    p = subprocess.run(["samtools","coverage","-Q",str(min_mapq),str(bam)], check=True, capture_output=True, text=True)
    out = {}
    for line in p.stdout.splitlines():
        if not line or line.startswith("#"): continue
        f = line.split("\t")
        if len(f) >= 9:
            out[f[0]] = {"mapped_reads":float(f[3]),"covered_bases":float(f[4]),"breadth":float(f[5])/100.0,"mean_depth":float(f[6]),"mean_base_quality":float(f[7]),"mean_mapping_quality":float(f[8])}
    return out

def low_depth_bed(bam: Path, path: Path, min_depth: float, min_mapq: int) -> None:
    p = subprocess.run(["samtools","depth","-aa","-d","0","-Q",str(min_mapq),str(bam)], check=True, capture_output=True, text=True)
    intervals=[]; cur=None
    for line in p.stdout.splitlines():
        chrom, pos, dep = line.split("\t")[:3]; pos0=int(pos)-1
        if float(dep) >= min_depth:
            if cur: intervals.append(cur); cur=None
            continue
        if cur and cur[0]==chrom and cur[2]==pos0: cur=(chrom,cur[1],pos0+1)
        else:
            if cur: intervals.append(cur)
            cur=(chrom,pos0,pos0+1)
    if cur: intervals.append(cur)
    with path.open("w") as h:
        for x in intervals: h.write(f"{x[0]}\t{x[1]}\t{x[2]}\n")

def consensus(reference: Path, bam: Path, prefix: Path, min_depth: float, min_mapq: int, basequal: int) -> Path:
    vcf=prefix.with_suffix(".vcf.gz"); mask=prefix.with_suffix(".low_depth.bed"); fa=prefix.with_suffix(".consensus.fasta")
    mp = subprocess.Popen(["bcftools","mpileup","-Ou","-f",str(reference),"-q",str(min_mapq),"-Q",str(basequal),"-d","100000",str(bam)], stdout=subprocess.PIPE)
    assert mp.stdout is not None
    call = subprocess.run(["bcftools","call","-mv","--ploidy","1","-Oz","-o",str(vcf)], stdin=mp.stdout)
    mp.stdout.close(); status=mp.wait()
    if status or call.returncode: raise subprocess.CalledProcessError(status or call.returncode, "bcftools")
    run(["bcftools","index","--force",str(vcf)]); low_depth_bed(bam, mask, min_depth, min_mapq)
    cmd=["bcftools","consensus","-f",str(reference)]
    if mask.stat().st_size: cmd += ["-m",str(mask)]
    cmd.append(str(vcf))
    with fa.open("w") as out: subprocess.run(cmd, check=True, stdout=out)
    return fa

def align_identity(reference: Path, consensus_fa: Path) -> dict[str, dict[str, float]]:
    p=subprocess.run(["minimap2","-x","asm5","--secondary=no","-c",str(reference),str(consensus_fa)], check=True, capture_output=True, text=True)
    best={}
    for line in p.stdout.splitlines():
        f=line.split("\t")
        if len(f)<12: continue
        q=f[0]; target=f[5]
        if q != target: continue
        matches=int(f[9]); block=int(f[10]); mapq=int(f[11]); cand={"identity":matches/block if block else None,"identical_positions":matches,"aligned_positions":block,"consensus_mapq":mapq,"identity_method":"minimap2_asm5"}
        if target not in best or cand["identical_positions"]>best[target]["identical_positions"]: best[target]=cand
    return best

def fixed_coordinate_identity(ref: str, cons: str) -> dict[str, object] | None:
    ref=ref.upper(); cons=cons.upper()
    n=min(len(ref),len(cons)); compared=0; same=0
    for a,b in zip(ref[:n],cons[:n]):
        if a=="N" or b=="N": continue
        compared += 1
        same += int(a==b)
    extra=abs(len(ref)-len(cons))
    compared += extra
    if compared == 0: return None
    return {"identity":same/compared,"identical_positions":same,"aligned_positions":compared,"identity_method":"fixed_coordinate_global"}

def classify_gene_evidence(*, mapped_reads: float, breadth: float, mean_depth: float, identity: float | None, min_breadth: float, min_depth: float, min_identity: float) -> dict[str, object]:
    if mapped_reads == 0 or breadth == 0:
        return {"validation_state":"not_detected","validated_call":0,"final_call_source":"read_validation","decision_reason":"mapped_reads=0 or breadth=0"}
    if mean_depth < min_depth:
        return {"validation_state":"low_depth","validated_call":"","final_call_source":"initial_call_unresolved","decision_reason":"breadth>0 and mean_depth below threshold"}
    if breadth < min_breadth:
        return {"validation_state":"partial_coverage","validated_call":0,"final_call_source":"read_validation","decision_reason":"depth passes but breadth below threshold"}
    if identity is None:
        return {"validation_state":"identity_unresolved","validated_call":"","final_call_source":"initial_call_unresolved","decision_reason":"breadth/depth pass but identity could not be measured"}
    if identity < min_identity:
        return {"validation_state":"divergent","validated_call":0,"final_call_source":"read_validation","decision_reason":"breadth/depth pass but identity below threshold"}
    return {"validation_state":"confirmed_present","validated_call":1,"final_call_source":"read_validation","decision_reason":"breadth/depth/identity pass"}

def validation_decision_logic_rows(min_breadth: str = "MIN_BREADTH", min_depth: str = "MIN_DEPTH", min_identity: str = "MIN_IDENTITY") -> list[list[str]]:
    return [
        ["not_detected","mapped_reads=0 OR breadth=0","0","No read evidence for the gene"],
        ["low_depth",f"breadth>0 AND mean_depth < {min_depth}","preserve initial call","Insufficient depth for a confident read-validation decision"],
        ["partial_coverage",f"mean_depth >= {min_depth} AND breadth < {min_breadth}","0","Reads cover only part of the gene; not an intact-gene call"],
        ["identity_unresolved","breadth/depth pass but identity could not be measured","preserve initial call","Coverage supports the locus but sequence identity is unresolved"],
        ["divergent",f"breadth/depth pass AND identity < {min_identity}","0","Covered sequence is below identity threshold"],
        ["confirmed_present",f"breadth >= {min_breadth} AND depth >= {min_depth} AND identity >= {min_identity}","1","Read evidence confirms the gene"],
        ["not_tested_carried_forward","gene not selected for validation","preserve initial call","Initial pangenome call was carried forward"],
    ]

def validate_isolate(reference: Path, key_tsv: Path, r1: str, r2: str, outdir: Path, threads: int, min_breadth: float, min_depth: float, min_identity: float, min_mapq: int, basequal: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True); bam=outdir/"gene_reads.bam"
    map_reads(reference,r1,r2,bam,threads,min_mapq,outdir/"bwa_mem.log")
    cov=coverage(bam,min_mapq); fa=consensus(reference,bam,outdir/"gene_reads",min_depth,min_mapq,basequal); aln=align_identity(reference,fa)
    refs=read_fasta(reference); cons=read_fasta(fa)
    with key_tsv.open(newline="") as h: keys=list(csv.DictReader(h,delimiter="\t"))
    rows=[]
    for row in keys:
        key=row["reference_id"]; gene=row["Gene"]; c=cov.get(key,{}); a=aln.get(key)
        breadth=float(c.get("breadth",0)); depth=float(c.get("mean_depth",0)); mapped=float(c.get("mapped_reads",0))
        if not a and mapped>0 and breadth>0:
            a=fixed_coordinate_identity(refs.get(key,""),cons.get(key,""))
        identity=None if not a else a.get("identity")
        decision=classify_gene_evidence(mapped_reads=mapped,breadth=breadth,mean_depth=depth,identity=identity,min_breadth=min_breadth,min_depth=min_depth,min_identity=min_identity)
        rows.append({"reference_id":key,"Gene":gene,**decision,"breadth":breadth,"percent_coverage":breadth*100.0,"mean_depth":depth,"identity":"NA" if identity is None else identity,"percent_identity":"NA" if identity is None else identity*100.0,"identity_method":a.get("identity_method","unresolved") if a else "unresolved","identical_positions":"" if not a else int(a.get("identical_positions",0)),"aligned_positions":"" if not a else int(a.get("aligned_positions",0)),"mapped_reads":int(mapped),"mean_mapping_quality":float(c.get("mean_mapping_quality",0))})
    fields=["reference_id","Gene","validated_call","validation_state","decision_reason","final_call_source","breadth","percent_coverage","mean_depth","identity","percent_identity","identity_method","identical_positions","aligned_positions","mapped_reads","mean_mapping_quality"]
    write_tsv(outdir/"metrics.tsv",fields,rows)
