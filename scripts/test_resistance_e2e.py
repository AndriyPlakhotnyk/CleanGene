#!/usr/bin/env python3
"""Real-tool execution test: paired reads -> Shovill -> Prokka -> Panaroo ->
read validation -> AMRFinderPlus -> operon reports, with CheckM2/downsampling off.
Synthetic in-silico fixtures test execution, not biological performance.
"""
from __future__ import annotations
import argparse
import gzip
import json
import random
import subprocess
import sys
from pathlib import Path

from cleangene.fasta import read_fasta
from cleangene.resistance import revcomp
from cleangene.util import read_tsv, sha256, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--amrfinder-db", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path, help="annotatable bacterial backbone, at least 200 kb")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    # Select a public tet(M) reference from the pinned AMRFinder database.
    refs = read_fasta(args.amrfinder_db / "AMR_CDS.fa")
    reference_id, reference = next((name, seq) for name, seq in refs.items() if "|tet(M)|" in name and len(seq) > 1800)
    backbone = max(read_fasta(args.genome).values(), key=len)[:200000]
    if len(backbone) < 200000:
        raise SystemExit("Backbone must contain at least 200 kb")
    rows = []
    for number, isolate in enumerate(("sample_a", "sample_b")):
        allele = reference
        if number:
            pos = 600
            allele = allele[:pos] + ("A" if allele[pos] != "A" else "C") + allele[pos+1:]
            allele = allele[:900] + "AAA" + allele[900:]
        sequence = backbone[:100000] + "TAA" * 20 + allele + "TAA" * 20 + backbone[100000:]
        paths = [root / f"{isolate}_R{n}.fastq.gz" for n in (1, 2)]
        rng = random.Random(300 + number)
        with gzip.open(paths[0], "wt") as one, gzip.open(paths[1], "wt") as two:
            for i in range(len(sequence)*200//300):
                start = rng.randrange(len(sequence)-400)
                for handle, read in ((one, sequence[start:start+150]), (two, revcomp(sequence[start+250:start+400]))):
                    handle.write(f"@r{i}\n{read}\n+\n{'I'*len(read)}\n")
        rows.append([isolate, "Enterococcus faecium", *map(str, paths)])
    manifest = root / "manifest.tsv"
    write_tsv(manifest, ["isolate_id", "group_id", "R1", "R2"], rows)
    cfg = {"TAXONOMY_MODE": "off", "AMRFINDER_DB": str(args.amrfinder_db.resolve()),
        "CPUS": str(args.threads), "VALIDATION_CPUS": str(args.threads), "ARBITRATION_CPUS": str(args.threads),
        "PANAROO_SMALL_CPUS": str(args.threads), "RESISTANCE_CPUS": str(args.threads), "RESISTANCE_MERGE_CPUS": str(args.threads),
        "SHOVILL_MEMORY_GB": "8", "PREPROCESS_USE_NODE_LOCAL_SCRATCH": "false"}
    config = root / "config.env"
    config.write_text("".join(f'{k}="{v}"\n' for k, v in cfg.items()))
    command = [sys.executable, "-m", "cleangene", "run", "--profile", "local", "--manifest", str(manifest),
        "--config", str(config), "--analysis-root", str(root), "--run-id", "e2e", "--skip-downsampling", "--ignore-checkm2", "-ressitanec-operon"]
    subprocess.run(command, check=True)
    run = root / "runs" / "e2e"
    resolved = json.loads((run / "provenance" / "resolved_config.json").read_text())
    assert resolved["SKIP_DOWNSAMPLING"] == "true" and resolved["CHECKM2_MODE"] == "off"
    qc = read_tsv(run / "results" / "cohort" / "isolate_qc.tsv")
    assert len(qc) == 2 and all(r["excluded"] == "0" for r in qc), qc
    assert all(not r.get("checkm2_completeness") for r in qc)
    for log in run.rglob("shovill.stderr"):
        text = log.read_text()
        assert "No read depth reduction requested or necessary" in text
        assert "seqkit sample" not in text
    out = run / "results" / "resistance_analysis"
    loci = read_tsv(out / "tables" / "loci.tsv")
    target = [l for l in loci if l["operon_type"] == "tet(M)"]
    assert len(target) == 2, loci
    assert len({l["variant_id"] for l in target}) == 2, target
    assert len({l["bin_id"] for l in target}) == 1, target
    assert all(float(l["read_supported_fraction"]) >= .95 for l in target), target
    assert all(l["origin"] == "unknown" for l in target)
    events = [e for e in read_tsv(out / "tables" / "variant_events.tsv") if e["operon_type"] == "tet(M)"]
    assert any(e["event"] == "SNP" for e in events), events
    assert any(e["event"] in {"insertion", "deletion"} for e in events), events
    assert list((out / "figures").rglob("alignment_sites.pdf"))
    assert read_tsv(run / "results" / "groups" / "Enterococcus_faecium" / "cleaned_pangenome.tsv")
    markers = {str(p): p.stat().st_mtime_ns for p in (run / "state" / "resistance_scan").glob("*.done.json")}
    before = sha256(out / "tables" / "loci.tsv")
    report_mtime = (out / "tables" / "loci.tsv").stat().st_mtime_ns
    subprocess.run([sys.executable, "-m", "cleangene", "run", "--profile", "local", "--analysis-root", str(root), "--resume", "e2e"], check=True)
    assert sha256(out / "tables" / "loci.tsv") == before
    assert (out / "tables" / "loci.tsv").stat().st_mtime_ns == report_mtime, "Resume regenerated unchanged reports"
    assert markers == {str(p): p.stat().st_mtime_ns for p in (run / "state" / "resistance_scan").glob("*.done.json")}
    report = {"status": "PASS", "run_dir": str(run), "fixture": "two in-silico 200-kb bacterial backbones with public tet(M) reference; not a biological E. faecium benchmark",
        "reference_id": reference_id, "isolates": 2, "resistance_loci": len(loci), "target_exact_variants": 2,
        "target_similarity_bins": 1, "SNP_and_indel_detected": True, "skip_downsampling": True, "checkm2": "off", "resume_cache": "PASS"}
    (root / "e2e_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
