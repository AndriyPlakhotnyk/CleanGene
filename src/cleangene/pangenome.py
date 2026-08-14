from __future__ import annotations
import csv, re
from pathlib import Path
from .fasta import read_fasta, write_fasta
from .util import safe_name, write_tsv

META = {"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"}

def present(value: str) -> int:
    return 0 if value.strip().lower() in {"", "0", "0.0", "na", "nan", "none", "-", "."} else 1

def normalize_panaroo(path: Path, isolates: list[str]) -> list[dict[str, object]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        column_for = {x: (x if x in fields else safe_name(x) if safe_name(x) in fields else "") for x in isolates}
        missing = [x for x, c in column_for.items() if not c]
        if missing: raise SystemExit("Panaroo matrix missing isolates: " + ", ".join(missing[:10]))
        rows = []
        seen: dict[str, int] = {}
        for line, row in enumerate(reader, 2):
            base = (row.get("Gene") or f"gene_row_{line}").strip()
            seen[base] = seen.get(base, 0) + 1
            gene = base if seen[base] == 1 else f"{base}__row{line}"
            rows.append({"Gene": gene, **{i: present(row.get(column_for[i], "")) for i in isolates}})
        return rows

def select_rows(rows: list[dict[str, object]], isolates: list[str], scope: str, cutoff: float) -> list[dict[str, object]]:
    if not 0 < cutoff <= 1: raise ValueError("accessory cutoff must be >0 and <=1")
    result = []
    for row in rows:
        n = sum(int(row[i]) for i in isolates); prevalence = n / len(isolates)
        if scope == "all" or (scope == "accessory" and n > 0 and prevalence <= cutoff) or (scope == "differential" and 0 < n < len(isolates)):
            result.append(row)
    return result

def recover_sequences(selected: list[dict[str, object]], panaroo_dir: Path) -> tuple[list[tuple[str, str]], list[list[object]]]:
    ref = read_fasta(panaroo_dir / "pan_genome_reference.fa")
    gene_data: dict[str, str] = {}
    gd = panaroo_dir / "gene_data.csv"
    if gd.is_file():
        with gd.open(newline="", errors="replace") as handle:
            for row in csv.DictReader(handle):
                locus = (row.get("annotation_id") or "").strip(); seq = (row.get("dna_sequence") or "").strip().upper()
                if locus and seq: gene_data[locus] = seq
    cluster_loci: dict[str, list[str]] = {}
    with (panaroo_dir / "gene_presence_absence.csv").open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            gene = (row.get("Gene") or "").strip(); loci = []
            for field, value in row.items():
                if field in META or not value: continue
                loci.extend(x for x in re.split(r"[;\s]+", value.strip()) if x)
            cluster_loci[gene] = loci
    records = []; sources = []
    for idx, row in enumerate(selected, 1):
        gene = str(row["Gene"]); base = gene.split("__row", 1)[0]
        seq = ref.get(gene) or ref.get(base); source = "pan_genome_reference"; locus = base
        if not seq:
            locus = next((x for x in cluster_loci.get(base, []) if x in gene_data), "")
            seq = gene_data.get(locus, ""); source = "panaroo_gene_data_member"
        if not seq: raise SystemExit(f"No recoverable sequence for Panaroo cluster {gene}")
        key = f"CG{idx:08d}"; records.append((key, seq)); sources.append([key, gene, source, locus, len(seq)])
    return records, sources

def write_binary(path: Path, rows: list[dict[str, object]], isolates: list[str]) -> None:
    write_tsv(path, ["Gene", *isolates], ([r["Gene"], *[int(r[i]) for i in isolates]] for r in rows))
