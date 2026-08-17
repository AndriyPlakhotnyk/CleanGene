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
        q=f[0]; matches=int(f[9]); block=int(f[10]); mapq=int(f[11]); cand={"identity":matches/block if block else 0.0,"identical_positions":matches,"aligned_positions":block,"consensus_mapq":mapq}
        if q not in best or cand["identical_positions"]>best[q]["identical_positions"]: best[q]=cand
    return best

def validate_isolate(reference: Path, key_tsv: Path, r1: str, r2: str, outdir: Path, threads: int, min_breadth: float, min_depth: float, min_identity: float, min_mapq: int, basequal: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True); bam=outdir/"gene_reads.bam"
    map_reads(reference,r1,r2,bam,threads,min_mapq,outdir/"bwa_mem.log")
    cov=coverage(bam,min_mapq); fa=consensus(reference,bam,outdir/"gene_reads",min_depth,min_mapq,basequal); aln=align_identity(reference,fa)
    with key_tsv.open(newline="") as h: keys=list(csv.DictReader(h,delimiter="\t"))
    rows=[]
    for row in keys:
        key=row["reference_id"]; gene=row["Gene"]; c=cov.get(key,{}); a=aln.get(key,{})
        breadth=float(c.get("breadth",0)); depth=float(c.get("mean_depth",0)); identity=float(a.get("identity",0)); call=int(breadth>=min_breadth and depth>=min_depth and identity>=min_identity)
        state="confirmed_present" if call else "partial_or_divergent" if breadth>0 else "not_detected"
        rows.append({"reference_id":key,"Gene":gene,"validated_call":call,"validation_state":state,"breadth":breadth,"percent_coverage":breadth*100.0,"mean_depth":depth,"identity":identity,"percent_identity":identity*100.0,"identical_positions":int(a.get("identical_positions",0)),"aligned_positions":int(a.get("aligned_positions",0)),"mapped_reads":int(c.get("mapped_reads",0)),"mean_mapping_quality":float(c.get("mean_mapping_quality",0))})
    fields=["reference_id","Gene","validated_call","validation_state","breadth","percent_coverage","mean_depth","identity","percent_identity","identical_positions","aligned_positions","mapped_reads","mean_mapping_quality"]
    write_tsv(outdir/"metrics.tsv",fields,rows)
