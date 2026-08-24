from __future__ import annotations
import gzip
from pathlib import Path

def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}; name = ""
    handle_cm = gzip.open(path, "rt", errors="replace") if path.suffix == ".gz" else path.open(errors="replace")
    with handle_cm as handle:
        for raw in handle:
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                name = line[1:].split()[0]; records.setdefault(name, [])
            elif name: records[name].append(line.upper())
    return {k: "".join(v) for k, v in records.items()}

def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80): handle.write(seq[i:i+80] + "\n")

def assembly_metrics(path: Path) -> dict[str, object]:
    seqs = list(read_fasta(path).values())
    lengths = sorted((len(x) for x in seqs), reverse=True)
    total = sum(lengths); half = total / 2; acc = 0; n50 = 0; l50 = 0
    for i, length in enumerate(lengths, 1):
        acc += length
        if acc >= half: n50, l50 = length, i; break
    joined = "".join(seqs)
    gc = sum(b in "GC" for b in joined); atgc = sum(b in "ACGT" for b in joined)
    return {"assembly_length": total, "contigs": len(lengths), "n50": n50, "l50": l50,
            "ambiguous_bases": sum(b not in "ACGT" for b in joined), "gc_fraction": gc / atgc if atgc else 0.0}
