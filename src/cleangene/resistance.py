"""Assembly-derived resistance loci, read evidence, and reproducible variant bins.

A candidate locus is not evidence of co-transcription. Unknown contig origin and
unobserved flanks are retained explicitly, never converted to negative calls.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from .config import truthy
from .fasta import read_fasta, write_fasta
from .util import atomic_json, load_json, read_tsv, run, safe_name, write_tsv

VERSION = "3"
DEFAULT_DEFINITIONS = {
    "vanA": ["vanR-A", "vanS-A", "vanH-A", "vanA", "vanX-A", "vanY-A", "vanZ-A"],
    "vanB": ["vanR-B", "vanS-B", "vanY-B", "vanW-B", "vanH-B", "vanB", "vanX-B"],
}
RC = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


def revcomp(sequence: str) -> str:
    return sequence.translate(RC)[::-1]


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:16]


def root_dir(run_dir: Path) -> Path:
    return run_dir / "results" / "resistance_analysis"


def definitions(cfg: dict) -> dict[str, list[str]]:
    result = dict(DEFAULT_DEFINITIONS)
    if cfg.get("RESISTANCE_OPERON_DEFINITIONS"):
        result.update(load_json(Path(cfg["RESISTANCE_OPERON_DEFINITIONS"])))
    assigned = {}
    for name, genes in result.items():
        if not isinstance(name, str) or not isinstance(genes, list) or not genes or any(not isinstance(g, str) or not g for g in genes):
            raise ValueError("Operon definitions must map type names to nonempty lists of gene symbols")
        for gene in genes:
            if gene in assigned and assigned[gene] != name:
                raise ValueError(f"Operon gene {gene} occurs in multiple definitions")
            assigned[gene] = name
    return result


def validate_config(cfg: dict) -> None:
    if not 0 < float(cfg["RESISTANCE_IDENTITY"]) <= 1:
        raise ValueError("RESISTANCE_IDENTITY must be in (0, 1]")
    for key in ("RESISTANCE_CPUS", "RESISTANCE_MAX_GAP", "RESISTANCE_TOP", "RESISTANCE_MIN_DEPTH"):
        if int(cfg[key]) < 1:
            raise ValueError(f"{key} must be positive")
    for key in ("RESISTANCE_MIN_BREADTH", "RESISTANCE_MIN_ALLELE_FRACTION"):
        if not 0 < float(cfg[key]) <= 1:
            raise ValueError(f"{key} must be in (0, 1]")
    if int(cfg["RESISTANCE_MIN_MAPQ"]) < 0 or int(cfg["BASEQUAL"]) < 0:
        raise ValueError("Mapping and base quality thresholds must be nonnegative")
    definitions(cfg)
    if cfg.get("RESISTANCE_CONTIG_ORIGINS"):
        origins(Path(cfg["RESISTANCE_CONTIG_ORIGINS"]))


def preflight(cfg: dict, out: Path) -> None:
    validate_config(cfg)
    required = ["amrfinder", "mafft"]
    if truthy(cfg["RESISTANCE_READ_SUPPORT"]):
        required += ["bwa", "samtools"]
    missing = [tool for tool in required if not shutil.which(tool)]
    if missing:
        raise RuntimeError("Resistance analysis requires: " + ", ".join(missing) + "; update environment.yml and provision an AMRFinderPlus database")
    import edlib  # noqa: F401
    import matplotlib  # noqa: F401
    out.mkdir(parents=True, exist_ok=True)
    command = ["amrfinder", "--database_version"]
    if cfg.get("AMRFINDER_DB"):
        command += ["--database", cfg["AMRFINDER_DB"]]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    (out / "amrfinder_version.txt").write_text("\n".join(line for line in (result.stdout + result.stderr).splitlines() if not line.startswith("amrfinder took")) + "\n")
    atomic_json(out / "method.json", {"version": VERSION, "definitions": definitions(cfg),
        "settings": {k: v for k, v in cfg.items() if k.startswith(("RESISTANCE_", "AMRFINDER_"))},
        "identity": "1 - global Levenshtein edit distance / max(ungapped sequence lengths); gaps and ambiguous bases penalized",
        "clustering": "deterministic abundance-ordered greedy representative; not transitive single-linkage",
        "variant_source": "assembly; read evidence measures support for assembly bases, not transcription or complete molecular phasing"})


def origins(path: Path) -> dict[tuple[str, str], dict]:
    result = {}
    for row in read_tsv(path):
        if not all(row.get(key) for key in ("isolate_id", "contig", "origin", "evidence")):
            raise ValueError("Contig origins TSV requires isolate_id, contig, origin, evidence")
        if row["origin"] not in {"chromosome", "plasmid", "unknown"}:
            raise ValueError("Contig origin must be chromosome, plasmid, or unknown")
        key = row["isolate_id"], row["contig"]
        if key in result:
            raise ValueError(f"Duplicate contig origin: {key}")
        result[key] = row
    return result


def features_from_gff(path: Path, assembly: dict[str, str]) -> list[dict]:
    from urllib.parse import unquote
    result = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            contig, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            if contig not in assembly or not 1 <= start <= end <= len(assembly[contig]):
                raise ValueError(f"GFF coordinates do not match assembly: {contig}:{start}-{end}")
            attrs = {k: unquote(v) for k, v in (part.split("=", 1) for part in fields[8].split(";") if "=" in part)}
            sequence = assembly[contig][start-1:end]
            result.append({"contig": contig, "start": start, "end": end, "strand": strand,
                "id": attrs.get("ID", attrs.get("locus_tag", f"{contig}:{start}-{end}")),
                "gene": attrs.get("gene", ""), "product": attrs.get("product", ""),
                "sequence": revcomp(sequence) if strand == "-" else sequence})
    if not result:
        raise ValueError(f"No CDS features in {path}; operon analysis requires annotation")
    return sorted(result, key=lambda f: (f["contig"], f["start"], f["end"], f["id"]))


def parse_hits(path: Path) -> list[dict]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        symbol_field = "Element symbol" if "Element symbol" in fields else "Gene symbol"
        if not {symbol_field, "Contig id", "Start", "Stop", "Strand", "Type", "Method"}.issubset(fields):
            raise ValueError(f"Unrecognized AMRFinderPlus output header: {fields}")
        hits = []
        for row in reader:
            if row["Type"] != "AMR":
                continue
            if row["Contig id"] in {"NA", "", "."}:
                raise ValueError("AMRFinderPlus hit has no nucleotide coordinates")
            start, end = sorted((int(row["Start"]), int(row["Stop"])))
            hits.append({"symbol": row[symbol_field], "contig": row["Contig id"], "start": start,
                "end": end, "strand": row["Strand"], "method": row["Method"],
                "subtype": row.get("Subtype", ""), "class": row.get("Class", ""),
                "protein": row.get("Protein id", ""), "raw": row})
    return hits


def gene_symbol(hit: dict) -> str:
    symbol = hit["symbol"]
    if hit.get("subtype", "").startswith("POINT") or hit.get("method", "").startswith("POINT"):
        # AMRFinder reports gene_mutation; content categories describe genes, not their alleles.
        return symbol.rsplit("_", 1)[0]
    return symbol


def extract_loci(isolate: str, group: str, assembly: dict[str, str], features: list[dict],
                 hits: list[dict], defs: dict[str, list[str]], max_gap: int) -> list[dict]:
    gene_type = {gene: name for name, genes in defs.items() for gene in genes}
    by_contig = defaultdict(list)
    for feature in features:
        by_contig[feature["contig"]].append(feature)
    grouped = defaultdict(list)
    for hit in hits:
        if hit["contig"] not in assembly or not 1 <= hit["start"] <= hit["end"] <= len(assembly[hit["contig"]]):
            raise ValueError(f"AMRFinder coordinates do not match assembly: {hit}")
        # Known point mutations are single-gene contexts, not acquired operons.
        gene = gene_symbol(hit)
        kind = gene_type.get(gene, gene)
        grouped[kind, hit["contig"]].append(hit)
    result = []
    for (kind, contig), seeds in sorted(grouped.items()):
        # Annotation may recover accessory operon members that AMRFinder does not report.
        annotated = []
        for feature in by_contig[contig]:
            if feature["gene"] in defs.get(kind, []):
                annotated.append({**feature, "symbol": feature["gene"], "method": "annotation"})
        unique = {(h["symbol"], h["start"], h["end"]): h for h in annotated + seeds}
        members = sorted(unique.values(), key=lambda h: (h["start"], h["end"], h["symbol"]))
        blocks = []
        for hit in members:
            repeated = bool(blocks) and any(gene_symbol(h) == gene_symbol(hit) and h["end"] < hit["start"] for h in blocks[-1])
            if not blocks or repeated or hit["start"] - max(h["end"] for h in blocks[-1]) > max_gap:
                blocks.append([])
            blocks[-1].append(hit)
        for block in blocks:
            actual = [h for h in seeds if h in block]
            if not actual:
                continue
            left, right = min(h["start"] for h in block), max(h["end"] for h in block)
            core = [f for f in by_contig[contig] if f["end"] >= left and f["start"] <= right]
            if core:
                left = min(left, min(f["start"] for f in core))
                right = max(right, max(f["end"] for f in core))
            preceding = [f for f in by_contig[contig] if f["end"] < left]
            following = [f for f in by_contig[contig] if f["start"] > right]
            before = max(preceding, key=lambda f: f["end"]) if preceding else None
            after = min(following, key=lambda f: f["start"]) if following else None
            start = before["start"] if before else left
            end = after["end"] if after else right
            anchor = next((h for h in actual if h["symbol"] == kind), actual[0])
            strand = anchor["strand"]
            chosen = ([before] if before else []) + core + ([after] if after else [])
            annotated_features = []
            symbols_by_feature = defaultdict(set)
            for hit in block:
                candidates = [f for f in core if min(hit["end"], f["end"]) >= max(hit["start"], f["start"])]
                if candidates:
                    matching_id = hit.get("protein") or hit.get("id")
                    best = max(candidates, key=lambda f: (f["id"] == matching_id,
                        min(hit["end"], f["end"]) - max(hit["start"], f["start"]) + 1, f["id"]))
                    symbols_by_feature[best["id"]].add(gene_symbol(hit))
            for feature in chosen:
                f = dict(feature)
                # Prefer the reported protein ID, then maximum overlap; do not label neighboring overlapping CDS as AMR.
                f["amr_symbols"] = sorted(symbols_by_feature[f["id"]])
                role = "upstream" if feature is before else "downstream" if feature is after else "core"
                if strand == "-" and role != "core":
                    role = "downstream" if role == "upstream" else "upstream"
                f["role"] = role
                f["orientation"] = "+" if f["strand"] == strand else "-"
                annotated_features.append(f)
            if strand == "-":
                annotated_features.reverse()
            observed = sorted({gene_symbol(h) for h in block})
            missing = sorted(set(defs.get(kind, [])) - set(observed))
            sequence = assembly[contig][start-1:end]
            if strand == "-":
                sequence = revcomp(sequence)
            partial = any("PARTIAL" in h["method"] or "INTERNAL_STOP" in h["method"] for h in actual)
            locus_id = "L" + digest([group, isolate, kind, contig, start, end])
            result.append({"locus_id": locus_id, "isolate_id": isolate, "group_id": group, "operon_type": kind,
                "boundary_method": "defined_gene_span" if kind in defs else "single_determinant_context",
                "contig": contig, "start": start, "end": end, "core_start": left, "core_end": right,
                "strand": strand, "sequence": sequence, "features": annotated_features,
                "observed_genes": observed, "reported_amr_symbols": sorted({h["symbol"] for h in actual}), "unobserved_expected_genes": missing,
                "upstream_missing": after is None if strand == "-" else before is None,
                "downstream_missing": before is None if strand == "-" else after is None,
                "partial_amr_hit": partial,
                "context_status": "incomplete" if before is None or after is None or missing or partial else "complete",
                "amr_methods": sorted({h["method"] for h in actual}),
                "origin": "unknown", "origin_evidence": "no contig-level origin evidence supplied"})
    return result


def pileup_counts(bases: str) -> tuple[int, int, int]:
    """Count assembly-matching bases, observed bases, and indel-bearing reads."""
    match = total = indels = i = 0
    while i < len(bases):
        char = bases[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            i += 1
            j = i
            while i < len(bases) and bases[i].isdigit():
                i += 1
            if j == i:
                raise ValueError("Malformed mpileup indel")
            i += int(bases[j:i])
            indels += 1
            continue
        if char in ".,":
            match += 1
            total += 1
        elif char in "ACGTNacgtn*#":
            total += 1
        # Reference skips (< >) are not evidence for this assembly base.
        i += 1
    return match, total, indels


def read_support(loci: list[dict], assembly_path: Path, row: dict, out: Path, cfg: dict) -> None:
    for locus in loci:
        locus.update(evidence="assembly_only", read_supported_fraction="NA", mean_depth="NA",
                     discordant_positions="NA", read_note="read support disabled or reads unavailable")
    if not loci or not truthy(cfg["RESISTANCE_READ_SUPPORT"]) or not row.get("R1"):
        return
    from .evidence import map_reads
    # Use the entire isolate assembly so homologues compete for mappings.
    reference = out / "mapping_reference.fasta"
    write_fasta(reference, list(read_fasta(assembly_path).items()))
    run(["bwa", "index", str(reference)], stderr=out / "bwa_index.log")
    bam = out / "reads.bam"
    map_reads(reference, row["R1"], row.get("R2", ""), bam, int(cfg["RESISTANCE_CPUS"]),
              int(cfg["RESISTANCE_MIN_MAPQ"]), out / "mapping.log")
    bed = out / "loci.bed"
    bed.write_text("".join(f"{l['contig']}\t{l['start']-1}\t{l['end']}\n" for l in loci))
    run(["samtools", "faidx", str(reference)])
    command = ["samtools", "mpileup", "-aa", "-B", "-d", "1000000", "-q", cfg["RESISTANCE_MIN_MAPQ"],
               "-Q", cfg["BASEQUAL"], "-f", str(reference), "-l", str(bed), str(bam)]
    positions = {}
    with (out / "mpileup.log").open("w") as log, gzip.open(out / "base_support.tsv.gz", "wt") as output:
        output.write("contig\tposition\tassembly_base\tdepth\tmatching_bases\tindel_reads\tsupported\n")
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True)
        assert process.stdout is not None
        try:
            for line in process.stdout:
                fields = line.rstrip().split("\t")
                if len(fields) < 6:
                    raise ValueError("Malformed mpileup output")
                contig, pos, base = fields[:3]
                matches, depth, indels = pileup_counts(fields[4]) if int(fields[3]) else (0, 0, 0)
                fraction = max(0, matches - indels) / depth if depth else 0
                supported = depth >= int(cfg["RESISTANCE_MIN_DEPTH"]) and fraction >= float(cfg["RESISTANCE_MIN_ALLELE_FRACTION"]) and base.upper() in "ACGT"
                discordant = depth >= int(cfg["RESISTANCE_MIN_DEPTH"]) and fraction < float(cfg["RESISTANCE_MIN_ALLELE_FRACTION"])
                positions[contig, int(pos)] = (depth, supported, discordant)
                output.write(f"{contig}\t{pos}\t{base}\t{depth}\t{matches}\t{indels}\t{int(supported)}\n")
        finally:
            process.stdout.close()
            code = process.wait()
        if code:
            raise subprocess.CalledProcessError(code, command)
    for locus in loci:
        values = [positions.get((locus["contig"], p), (0, False, False)) for p in range(locus["start"], locus["end"] + 1)]
        fraction = sum(v[1] for v in values) / len(values)
        discordant = sum(v[2] for v in values)
        locus.update(read_supported_fraction=fraction, mean_depth=sum(v[0] for v in values) / len(values),
                     discordant_positions=discordant, evidence="assembly_and_reads" if fraction >= float(cfg["RESISTANCE_MIN_BREADTH"]) and not discordant else "assembly_read_warning",
                     read_note="assembly allele support; low-MAPQ repeats unresolved; short reads do not phase complete loci")
    # Keep the compact base-level evidence, avoiding an extra genome BAM per isolate.
    if not truthy(cfg["RESISTANCE_KEEP_BAM"]):
        for path in [bam, Path(str(bam) + ".bai"), bam.with_suffix(".mapping.json"), *out.glob("mapping_reference.fasta*")]:
            path.unlink(missing_ok=True)


def sample_inputs(run_dir: Path, index: int, cfg: dict) -> tuple[dict, bool, str]:
    from .workers import task_row, find_isolate_qc, user_excluded, manifest_row_for_task
    row = task_row(run_dir, "isolate", index)
    if "R1" not in row:
        row = manifest_row_for_task(row, read_tsv(run_dir / "provenance" / "manifest.tsv"))
    iso = row["isolate_id"]
    qc = read_tsv(find_isolate_qc(run_dir, row))[0]
    excluded = user_excluded(row) or qc.get("PASS/FAIL") == "FAIL" or qc["excluded"] in {"1", "true", "True"}
    row = {**row, **{k: qc[k] for k in ("assembly", "gff", "R1", "R2") if k in qc}}
    protein = Path(row.get("protein_fasta") or Path(row.get("gff", "")).with_suffix(".faa")) if row.get("gff") else Path("missing")
    input_paths = [Path(row.get(k, "")) for k in ("assembly", "gff", "R1", "R2")] + [protein, Path(str(protein) + ".gz")]
    signature = digest([VERSION, {k: v for k, v in cfg.items() if k.startswith(("RESISTANCE_", "AMRFINDER_"))},
        definitions(cfg), cfg["BASEQUAL"], excluded, row, [(str(p), p.stat().st_size, p.stat().st_mtime_ns)
        for p in input_paths if p.is_file()],
        (root_dir(run_dir) / "provenance" / "amrfinder_version.txt").read_text()])
    return row, excluded, signature


def isolate_task(run_dir: Path, index: int) -> None:
    from .workers import context
    cfg, _ = context(run_dir)
    row, excluded, signature = sample_inputs(run_dir, index, cfg)
    iso = row["isolate_id"]
    out = root_dir(run_dir) / "isolates" / safe_name(iso)
    out.mkdir(parents=True, exist_ok=True)
    marker = run_dir / "state" / "resistance_scan" / f"{safe_name(iso)}.done.json"
    if marker.is_file() and load_json(marker).get("signature") == signature and scan_outputs_present(out, excluded, cfg):
        return
    marker.unlink(missing_ok=True)
    if excluded:
        atomic_json(out / "loci.json", {"status": "excluded", "isolate_id": iso, "group_id": row["group_id"], "loci": []})
        atomic_json(marker, {"status": "complete", "signature": signature})
        return
    assembly_path, gff = Path(row.get("assembly", "")), Path(row.get("gff", ""))
    if not assembly_path.is_file() or not gff.is_file():
        raise RuntimeError(f"Resistance analysis requires an assembly and GFF for retained isolate {iso}")
    assembly = read_fasta(assembly_path)
    features = features_from_gff(gff, assembly)
    amr_assembly = out / "amr_assembly.fasta"
    write_fasta(amr_assembly, list(assembly.items()))
    command = ["amrfinder", "--nucleotide", str(amr_assembly), "--threads", cfg["RESISTANCE_CPUS"],
               "--plus", "--output", str(out / "amrfinder.tsv")]
    # Protein search improves sensitivity; nucleotide-only remains supported for compressed/non-Prokka annotations.
    protein = Path(row.get("protein_fasta") or gff.with_suffix(".faa"))
    compressed_protein = Path(str(protein) + ".gz") if protein.suffix != ".gz" else protein
    if compressed_protein.is_file() and (not protein.is_file() or protein.suffix == ".gz"):
        with gzip.open(compressed_protein, "rt") as source:
            protein = out / "amr_proteins.faa"
            protein.write_text(source.read())
    if protein.is_file() and cfg["AMRFINDER_ANNOTATION_FORMAT"] == "prokka":
        command += ["--protein", str(protein), "--gff", str(gff), "--annotation_format", "prokka"]
    organism = cfg.get("AMRFINDER_ORGANISM", "").strip()
    if not organism and row["group_id"].replace("_", " ").lower() == "enterococcus faecium":
        organism = "Enterococcus_faecium"
    if organism:
        command += ["--organism", organism]
    if cfg.get("AMRFINDER_DB"):
        command += ["--database", cfg["AMRFINDER_DB"]]
    atomic_json(out / "command.json", command)
    run(command, stderr=out / "amrfinder.log")
    amr_assembly.unlink()
    (out / "amr_proteins.faa").unlink(missing_ok=True)
    hits = parse_hits(out / "amrfinder.tsv")
    loci = extract_loci(iso, row["group_id"], assembly, features, hits, definitions(cfg), int(cfg["RESISTANCE_MAX_GAP"]))
    read_support(loci, assembly_path, row, out, cfg)
    atomic_json(out / "loci.json", {"status": "evaluated", "isolate_id": iso, "group_id": row["group_id"], "loci": loci})
    atomic_json(marker, {"status": "complete", "signature": signature, "loci": len(loci)})


def global_identity(a: str, b: str, minimum: float = 0) -> float:
    import edlib
    length = max(len(a), len(b))
    if not length:
        return 1.0
    if min(len(a), len(b)) / length < minimum:
        return 0.0
    # Distinct replacement symbols ensure an N/N pair cannot count as a match.
    a = re.sub("[^ACGT]", "X", a.upper())
    b = re.sub("[^ACGT]", "Y", b.upper())
    limit = math.floor((1 - minimum) * length + 1e-8) if minimum else -1
    distance = edlib.align(a, b, mode="NW", task="distance", k=limit)["editDistance"]
    return 0.0 if distance < 0 else 1 - distance / length


def cluster_sequences(sequences: dict[str, str], weights: dict[str, int], minimum: float) -> dict[str, str]:
    representatives = []
    assignments = {}
    for name in sorted(sequences, key=lambda n: (-weights.get(n, 0), -len(sequences[n]), n)):
        representative = next((r for r in representatives if global_identity(sequences[name], sequences[r], minimum) + 1e-12 >= minimum), None)
        if representative is None:
            representative = name
            representatives.append(name)
        assignments[name] = representative
    return assignments


def panaroo_families(run_dir: Path, group: str, isolates: list[str]) -> dict[tuple[str, str], str]:
    from .workers import prepared_pangenome_dir
    root = run_dir / "results" / "groups" / safe_name(group)
    path = prepared_pangenome_dir(run_dir, group, root) / "gene_presence_absence.csv"
    result = {}
    if path.is_file():
        with path.open() as handle:
            for row in csv.DictReader(handle):
                for iso in isolates:
                    for locus in re.split(r"[;\s]+", row.get(iso, row.get(safe_name(iso), "")) or ""):
                        if locus:
                            result[iso, locus] = "pan:" + row["Gene"]
    return result


def assign_families(loci: list[dict], family_map: dict) -> None:
    missing = {}
    for locus in loci:
        for feature in locus["features"]:
            family = ";".join("amr:" + s for s in feature["amr_symbols"])
            family = family or family_map.get((locus["isolate_id"], feature["id"]), "")
            feature["family"] = family
            if not family:
                sequence = feature["sequence"]
                missing[digest(sequence)] = sequence
    fallback = cluster_sequences(missing, {}, .95)
    for locus in loci:
        for feature in locus["features"]:
            if not feature["family"]:
                feature["family"] = "seq95:" + fallback[digest(feature["sequence"])]
        # Presence/absence is separate from copy number, gene order, orientation, and role.
        content = sorted({f["family"] for f in locus["features"]} | {"amr:" + g for g in locus["observed_genes"]})
        if locus["upstream_missing"]:
            content.append("UNKNOWN_UPSTREAM")
        if locus["downstream_missing"]:
            content.append("UNKNOWN_DOWNSTREAM")
        locus["gene_content"] = sorted(content)
        locus["category_id"] = "C" + digest([locus["operon_type"], sorted(content)])
        locus["gene_order"] = [f"{f['role']}:{f['family']}:{f['orientation']}" for f in locus["features"]]
        locus["variant_id"] = "V" + digest([locus["operon_type"], locus["category_id"], locus["sequence"]])


def align_sequences(sequences: dict[str, str], path: Path, threads: int) -> dict[str, str]:
    raw = path.with_name(path.stem + ".unaligned.fasta")
    write_fasta(raw, sorted(sequences.items()))
    if len(sequences) == 1:
        write_fasta(path, list(sequences.items()))
    else:
        # FFT-NS-1 avoids --auto choosing a costly strategy for long multi-gene loci.
        run(["mafft", "--nuc", "--thread", str(threads), "--retree", "1", "--maxiterate", "0", "--inputorder", str(raw)],
            stdout=path, stderr=path.with_suffix(".log"))
    aligned = read_fasta(path)
    if set(aligned) != set(sequences) or len({len(s) for s in aligned.values()}) != 1:
        raise RuntimeError("MAFFT output missing sequences or has inconsistent alignment lengths")
    if any(aligned[n].replace("-", "") != s for n, s in sequences.items()):
        raise RuntimeError("MAFFT changed input sequence content")
    return aligned


def variant_events(reference: str, alternate: str) -> list[dict]:
    if len(reference) != len(alternate):
        raise ValueError("Variant calls require equal aligned lengths")
    events = []
    refpos = altpos = 0
    current = None
    for column, (r, a) in enumerate(zip(reference, alternate), 1):
        rp, ap = refpos + (r != "-"), altpos + (a != "-")
        kind = "" if r == a else "uncertain" if r not in "ACGT-" or a not in "ACGT-" else "insertion" if r == "-" else "deletion" if a == "-" else "SNP"
        if kind:
            if current is None or current["event"] != kind or current["alignment_end"] != column - 1 or kind == "SNP":
                current = {"event": kind, "alignment_start": column, "alignment_end": column,
                           "reference_position": rp, "alternate_position": ap, "ref": "", "alt": ""}
                events.append(current)
            current["alignment_end"] = column
            current["ref"] += r if r != "-" else ""
            current["alt"] += a if a != "-" else ""
        else:
            current = None
        refpos, altpos = rp, ap
    for event in events:
        event["size"] = max(len(event["ref"]), len(event["alt"]))
        event["scale"] = "large" if event["size"] >= 50 else "small"
    return events


def plot_type(directory: Path, label: str, aligned: dict[str, str], ranked: list[str],
              counts: dict[str, int], variants: dict[str, dict], category_counts: dict[str, int],
              denominator: int, top: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    directory.mkdir(parents=True, exist_ok=True)
    chosen = ranked[:top]
    fig, axes = plt.subplots(1, 2, figsize=(15, max(4, len(chosen)*.38+1)))
    cats = sorted(category_counts, key=lambda c: (-category_counts[c], c))[:20]
    for ax, names, values, title in (
        (axes[0], cats, [category_counts[c] for c in cats], "Gene-content categories (top 20)"),
        (axes[1], chosen, [counts[v] for v in chosen], f"Exact variants (top {top})"),
    ):
        ax.barh(range(len(names)), [100*v/denominator for v in values], color="#28788e")
        ax.set_yticks(range(len(names)), names, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlim(0, max([100*v/denominator for v in values] + [1]) * 1.22)
        for i, value in enumerate(values):
            ax.text(100*value/denominator, i, f"  {value}/{denominator}", va="center", fontsize=8)
        ax.set_xlabel("% of evaluated isolates (copies counted once)")
        ax.set_title(title)
    fig.suptitle(label + " — prevalence; categories/variants can coexist in an isolate")
    fig.tight_layout()
    fig.savefig(directory / "prevalence.png", dpi=180)
    fig.savefig(directory / "prevalence.pdf")
    plt.close(fig)
    ref = aligned[ranked[0]]
    colors = ["#f2f2f2", "#e69f00", "#009e73", "#cc79a7", "#999999"]
    codes = {"match": 0, "SNP": 1, "insertion": 2, "deletion": 3, "uncertain": 4}
    def code(r, a):
        return codes["match" if r == a else "uncertain" if r not in "ACGT-" or a not in "ACGT-" else "insertion" if r == "-" else "deletion" if a == "-" else "SNP"]
    data = np.array([[code(r, a) for r, a in zip(ref, aligned[name])] for name in chosen])
    labels = [f"{v} | n={counts[v]} | {variants[v]['bin_id']}" for v in chosen]
    # A single-pixel SNP can disappear in a whole-locus raster: detailed pages below show every variable column.
    fig, ax = plt.subplots(figsize=(18, max(3, len(chosen)*.45+1.8)))
    ax.imshow(data, aspect="auto", interpolation="nearest", cmap=ListedColormap(colors), vmin=0, vmax=4)
    ax.set_yticks(range(len(chosen)), labels, fontsize=7)
    ax.set_xlabel("MSA column (gaps included); base-resolution differences in alignment_sites.pdf")
    ax.set_title(label + " — top exact variants vs most prevalent assembly haplotype")
    # Ensure isolated single-base events remain visible at overview resolution.
    for y, name in enumerate(chosen):
        differences = [i for i, value in enumerate(data[y]) if value]
        if len(differences) < 1000:
            ax.vlines(differences, y-.45, y+.45, colors=[colors[data[y, i]] for i in differences], linewidth=.7)
    fig.legend(handles=[Patch(color=c, label=n) for n, c in zip(codes, colors)], loc="lower center", ncol=5, fontsize=9)
    fig.tight_layout(rect=(0, .12, 1, 1))
    fig.savefig(directory / "alignment_overview.png", dpi=200)
    fig.savefig(directory / "alignment_overview.pdf")
    plt.close(fig)
    sites = [i for i in range(len(ref)) if any(aligned[v][i] != ref[i] for v in chosen)]
    reference_positions = []
    position = 0
    for char in ref:
        position += char != "-"
        reference_positions.append(str(position) + ("+" if char == "-" else ""))
    with PdfPages(directory / "alignment_sites.pdf") as pdf:
        for offset in range(0, max(1, len(sites)), 60):
            columns = sites[offset:offset+60]
            fig, ax = plt.subplots(figsize=(18, max(3, len(chosen)*.4+2)))
            if columns:
                ax.imshow(data[:, columns], aspect="auto", interpolation="nearest", cmap=ListedColormap(colors), vmin=0, vmax=4)
                for y, name in enumerate(chosen):
                    for x, col in enumerate(columns):
                        ax.text(x, y, aligned[name][col], ha="center", va="center", fontsize=7, family="monospace")
                ax.set_yticks(range(len(chosen)), labels, fontsize=6)
                ax.set_xticks(range(len(columns)), [f"{i+1}/{reference_positions[i]}" for i in columns], rotation=90, fontsize=6)
                ax.set_xlabel("MSA column / reference nucleotide (0+ is insertion before base 1; + is insertion after base)")
            else:
                ax.text(.5, .5, "No nucleotide differences among displayed variants", ha="center", va="center")
                ax.axis("off")
            ax.set_title(f"{label} — differing columns {offset+1}–{offset+len(columns)} of {len(sites)}; assembly-derived")
            fig.tight_layout()
            pdf.savefig(fig)
            if offset == 0:
                fig.savefig(directory / "alignment_sites_page1.png", dpi=180)
            plt.close(fig)
    fig, ax = plt.subplots(figsize=(18, max(3, len(chosen)*.8+1.5)))
    families = sorted({f["family"] for v in chosen for f in variants[v]["features"]})
    palette = plt.get_cmap("tab20")
    family_color = {f: palette(i % 20) for i, f in enumerate(families)}
    for y, name in enumerate(chosen):
        locus = variants[name]
        ax.plot([0, len(locus["sequence"])], [y, y], color="#666666", lw=1)
        for feature in locus["features"]:
            start = feature["start"]-locus["start"] if locus["strand"] == "+" else locus["end"]-feature["end"]
            length = feature["end"]-feature["start"]+1
            direction = 1 if feature["orientation"] == "+" else -1
            ax.arrow(start if direction == 1 else start+length, y, direction*length, 0,
                width=.14, head_width=.28, head_length=min(length*.2, 100), length_includes_head=True,
                color=family_color[feature["family"]], linewidth=.5)
            name_label = ",".join(feature["amr_symbols"]) or feature["gene"] or feature["family"]
            ax.text(start+length/2, y-.20, name_label, fontsize=7, ha="center", rotation=15)
    ax.set_yticks(range(len(chosen)), [f"{v} (n={counts[v]})" for v in chosen], fontsize=7)
    ax.set_ylim(len(chosen)-.5, -.7)
    ax.set_xlabel("Oriented assembly locus coordinate (includes one flanking CDS on each side)")
    ax.set_title(label + " — gene order, orientation and spacing; colours identify families")
    fig.tight_layout()
    fig.savefig(directory / "gene_structure.pdf")
    fig.savefig(directory / "gene_structure.png", dpi=180)
    plt.close(fig)


def merge(run_dir: Path) -> None:
    from .workers import context
    cfg, _ = context(run_dir)
    out = root_dir(run_dir)
    tables, figures, alignments = out / "tables", out / "figures", out / "alignments"
    marker = run_dir / "state" / "resistance_merge.done.json"
    tasks = read_tsv(run_dir / "state" / "isolate_tasks.tsv")
    from .workers import prepared_pangenome_dir
    inputs = []
    for task in tasks:
        inputs += [out / "isolates" / safe_name(task["isolate_id"]) / "loci.json",
                   run_dir / "state" / "resistance_scan" / f"{safe_name(task['isolate_id'])}.done.json"]
    for group in sorted({t["group_id"] for t in tasks}):
        inputs.append(prepared_pangenome_dir(run_dir, group, run_dir / "results/groups" / safe_name(group)) / "gene_presence_absence.csv")
    origin_data = origins(Path(cfg["RESISTANCE_CONTIG_ORIGINS"])) if cfg.get("RESISTANCE_CONTIG_ORIGINS") else {}
    signature = digest([VERSION, definitions(cfg), {k: v for k, v in cfg.items() if k.startswith(("RESISTANCE_", "AMRFINDER_"))},
        sorted((list(k), v) for k, v in origin_data.items()),
        [(str(p), p.stat().st_size, p.stat().st_mtime_ns) for p in inputs if p.is_file()]])
    if marker.is_file():
        previous = load_json(marker)
        if previous.get("signature") == signature and previous.get("outputs") and all((out / p).is_file() for p in previous["outputs"]):
            return
    marker.unlink(missing_ok=True)
    # These three directories contain only generated reports, and are replaced together.
    for path in (tables, figures, alignments):
        if path.is_dir():
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)
    samples, loci = [], []
    for task in tasks:
        path = out / "isolates" / safe_name(task["isolate_id"]) / "loci.json"
        done = run_dir / "state" / "resistance_scan" / f"{safe_name(task['isolate_id'])}.done.json"
        if not path.is_file() or not done.is_file():
            raise RuntimeError(f"Incomplete resistance task: {task['isolate_id']}")
        result = load_json(path)
        samples.append({k: result[k] for k in ("isolate_id", "group_id", "status")})
        loci.extend(result["loci"])
    supplied_origins = origin_data
    for locus in loci:
        hit = supplied_origins.get((locus["isolate_id"], locus["contig"]))
        if hit:
            locus["origin"], locus["origin_evidence"] = hit["origin"], hit["evidence"]
    write_tsv(tables / "sample_status.tsv", ["isolate_id", "group_id", "status"], samples)
    prevalence, variant_rows, category_rows, event_rows, gene_rows, structural_rows, bin_rows = [], [], [], [], [], [], []
    for group in sorted({s["group_id"] for s in samples}):
        evaluated = [s["isolate_id"] for s in samples if s["group_id"] == group and s["status"] == "evaluated"]
        cohort = [l for l in loci if l["group_id"] == group]
        assign_families(cohort, panaroo_families(run_dir, group, evaluated))
        for kind in sorted({l["operon_type"] for l in cohort}):
            members = [l for l in cohort if l["operon_type"] == kind]
            typeslug = safe_name(kind) + "_" + digest(kind)[:6]
            group_slug = safe_name(group)
            cats = defaultdict(list)
            for locus in members:
                cats[locus["category_id"]].append(locus)
            for cat, cat_members in sorted(cats.items()):
                sequences = {l["variant_id"]: l["sequence"] for l in cat_members}
                weights = {v: len({l["isolate_id"] for l in cat_members if l["variant_id"] == v}) for v in sequences}
                assignments = cluster_sequences(sequences, weights, float(cfg["RESISTANCE_IDENTITY"]))
                for locus in cat_members:
                    representative = assignments[locus["variant_id"]]
                    locus["bin_id"] = "B" + digest([cat, representative, cfg["RESISTANCE_IDENTITY"]])
                    locus["bin_representative"] = representative
                    locus["identity_to_representative"] = global_identity(locus["sequence"], sequences[representative])
                category_rows.append({"group_id": group, "operon_type": kind, "category_id": cat,
                    "gene_content": ";".join(cat_members[0]["gene_content"]),
                    "n_isolates": len({l["isolate_id"] for l in cat_members}), "n_copies": len(cat_members),
                    "n_bins": len(set(assignments.values()))})
                for representative in sorted(set(assignments.values())):
                    bin_members = [l for l in cat_members if l["bin_representative"] == representative]
                    bin_rows.append({"group_id": group, "operon_type": kind, "category_id": cat,
                        "bin_id": bin_members[0]["bin_id"], "representative_variant": representative,
                        "n_isolates": len({l["isolate_id"] for l in bin_members}), "n_copies": len(bin_members),
                        "n_exact_variants": len({l["variant_id"] for l in bin_members}),
                        "minimum_identity_to_representative": min(l["identity_to_representative"] for l in bin_members)})
            variants = {l["variant_id"]: l for l in sorted(members, key=lambda l: (l["isolate_id"], l["locus_id"]), reverse=True)}
            counts = {v: len({l["isolate_id"] for l in members if l["variant_id"] == v}) for v in variants}
            ranked = sorted(variants, key=lambda v: (-counts[v], v))
            target = alignments / group_slug / typeslug
            target.mkdir(parents=True, exist_ok=True)
            aligned = align_sequences({v: variants[v]["sequence"] for v in ranked}, target / "all_unique.aligned.fasta", int(cfg["RESISTANCE_MERGE_CPUS"]))
            ref_id = ranked[0]
            for rank, vid in enumerate(ranked, 1):
                locus = variants[vid]
                copies = [l for l in members if l["variant_id"] == vid]
                events = variant_events(aligned[ref_id], aligned[vid])
                event_counts = Counter(e["event"] for e in events)
                variant_rows.append({"group_id": group, "operon_type": kind, "rank": rank, "variant_id": vid,
                    "category_id": locus["category_id"], "bin_id": locus["bin_id"],
                    "n_isolates": counts[vid], "n_copies": len(copies), "n_evaluated": len(evaluated),
                    "prevalence": counts[vid]/len(evaluated), "length": len(locus["sequence"]),
                    "reference_variant": ref_id, "SNPs": event_counts["SNP"], "insertions": event_counts["insertion"],
                    "deletions": event_counts["deletion"], "uncertain_events": event_counts["uncertain"],
                    "read_supported_isolates": len({l["isolate_id"] for l in copies if l["evidence"] == "assembly_and_reads"}),
                    "origins": ";".join(sorted({l["origin"] for l in copies})),
                    "context_statuses": ";".join(sorted({l["context_status"] for l in copies}))})
                for event in events:
                    event_rows.append({"group_id": group, "operon_type": kind, "variant_id": vid,
                        "reference_variant": ref_id, "evidence": "assembly_alignment", **event})
                reference_locus = variants[ref_id]
                structural_rows.append({"group_id": group, "operon_type": kind, "variant_id": vid,
                    "reference_variant": ref_id, "gene_order": ";".join(locus["gene_order"]),
                    "gene_order_or_orientation_changed": int(locus["gene_order"] != reference_locus["gene_order"]),
                    "gene_copy_counts": json.dumps(dict(Counter(f["family"] for f in locus["features"])), sort_keys=True),
                    "length_difference": len(locus["sequence"]) - len(reference_locus["sequence"])})
            n_present = len({l["isolate_id"] for l in members})
            prevalence.append({"group_id": group, "operon_type": kind, "n_present": n_present,
                "n_evaluated": len(evaluated), "n_not_detected": len(evaluated)-n_present,
                "prevalence": n_present/len(evaluated), "n_copies": len(members),
                "n_categories": len(cats), "n_bins": len({l["bin_id"] for l in members}), "n_exact_variants": len(variants),
                "n_complete_context_isolates": len({l["isolate_id"] for l in members if l["context_status"] == "complete"})})
            category_counts = {c: len({l["isolate_id"] for l in rows}) for c, rows in cats.items()}
            plot_type(figures / group_slug / typeslug, f"{group}: {kind}", aligned, ranked, counts, variants,
                      category_counts, len(evaluated), int(cfg["RESISTANCE_TOP"]))
            # One binary matrix per type includes all inferred genes AND both flank families.
            genes = sorted({g for l in members for g in l["gene_content"]} | {"amr:" + g for g in definitions(cfg).get(kind, [])})
            write_tsv(tables / group_slug / typeslug / "gene_presence_absence.tsv", ["locus_id", "isolate_id", *genes],
                [[l["locus_id"], l["isolate_id"], *[int(g in l["gene_content"]) for g in genes]] for l in members])
    for locus in loci:
        for order, feature in enumerate(locus["features"], 1):
            gene_rows.append({"locus_id": locus["locus_id"], "isolate_id": locus["isolate_id"], "order": order,
                **{k: feature[k] for k in ("id", "family", "gene", "product", "role", "start", "end", "strand", "orientation")}})
    locus_fields = ["locus_id", "isolate_id", "group_id", "operon_type", "boundary_method", "category_id", "bin_id", "variant_id",
        "bin_representative", "identity_to_representative", "contig", "start", "end", "core_start", "core_end", "strand",
        "origin", "origin_evidence", "context_status", "partial_amr_hit", "upstream_missing", "downstream_missing", "observed_genes",
        "unobserved_expected_genes", "reported_amr_symbols", "amr_methods", "evidence", "read_supported_fraction", "mean_depth", "discordant_positions", "read_note"]
    write_tsv(tables / "loci.tsv", locus_fields,
              [{k: ";".join(map(str, l[k])) if isinstance(l.get(k), list) else l.get(k, "") for k in locus_fields} for l in loci])
    def table(name, fields, rows):
        write_tsv(tables / name, fields.split(), rows)
    table("operon_prevalence.tsv", "group_id operon_type n_present n_evaluated n_not_detected prevalence n_copies n_categories n_bins n_exact_variants n_complete_context_isolates", prevalence)
    table("categories.tsv", "group_id operon_type category_id gene_content n_isolates n_copies n_bins", category_rows)
    table("similarity_bins.tsv", "group_id operon_type category_id bin_id representative_variant n_isolates n_copies n_exact_variants minimum_identity_to_representative", bin_rows)
    table("variants.tsv", "group_id operon_type rank variant_id category_id bin_id n_isolates n_copies n_evaluated prevalence length reference_variant SNPs insertions deletions uncertain_events read_supported_isolates origins context_statuses", variant_rows)
    table("variant_events.tsv", "group_id operon_type variant_id reference_variant event alignment_start alignment_end reference_position alternate_position ref alt size scale evidence", event_rows)
    table("gene_structure.tsv", "group_id operon_type variant_id reference_variant gene_order gene_order_or_orientation_changed gene_copy_counts length_difference", structural_rows)
    table("locus_genes.tsv", "locus_id isolate_id order id family gene product role start end strand orientation", gene_rows)
    # Include zero-hit isolates in denominators; absence of detection is not proof of absence.
    type_columns = sorted({(l["group_id"], l["operon_type"]) for l in loci})
    detected = {(l["isolate_id"], l["group_id"], l["operon_type"]) for l in loci}
    write_tsv(tables / "isolate_presence_absence.tsv", ["isolate_id", "group_id", "status", *[f"{g}|{k}" for g, k in type_columns]],
        [[s["isolate_id"], s["group_id"], s["status"], *[int((s["isolate_id"], g, k) in detected) if s["group_id"] == g and s["status"] == "evaluated" else "NA" for g, k in type_columns]] for s in samples])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, max(3, len(prevalence)*.4+1)))
    if prevalence:
        ordered = sorted(prevalence, key=lambda r: (-r["prevalence"], r["group_id"], r["operon_type"]))
        ax.barh(range(len(ordered)), [100*r["prevalence"] for r in ordered], color="#28788e")
        ax.set_yticks(range(len(ordered)), [f"{r['group_id']}: {r['operon_type']} ({r['n_present']}/{r['n_evaluated']})" for r in ordered], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, 105)
        ax.set_xlabel("% of evaluated isolates with a detected resistance locus")
    else:
        ax.text(.5, .5, "No AMR loci detected in evaluated assemblies", ha="center", va="center")
        ax.axis("off")
    ax.set_title("Resistance locus prevalence (includes partial contexts)")
    fig.tight_layout()
    fig.savefig(figures / "operon_prevalence.png", dpi=180)
    fig.savefig(figures / "operon_prevalence.pdf")
    plt.close(fig)
    (out / "README.txt").write_text(
        "Assembly-derived candidate resistance operons/contexts; not transcriptional evidence.\n"
        "Gene-content categories include internal CDS and one upstream/downstream CDS. Missing flanks are UNKNOWN.\n"
        "95% bins use abundance-ordered global edit similarity to a fixed representative within each category.\n"
        "Exact sequence variants remain separate within bins; all are aligned, the top 12 per type are displayed by default.\n"
        "Events are relative to the most prevalent assembly sequence, not an ancestral or susceptible reference.\n"
        "Base-support TSVs retain read evidence. Assembly-only and read-warning variants remain visible.\n"
        "Small indels <50 bp; large >=50 bp. Gene order/copy differences are reported separately.\n"
        "MSA does not establish a rearrangement mechanism; repeat-rich/large changes require review.\n"
        "Short reads do not prove complete locus phasing. Split contigs are not joined.\n"
        "Unknown origin means chromosome/plasmid was not resolved; supplied origins retain evidence labels.\n"
        "Prevalence counts distinct retained/evaluated isolates; no-hit is not-detected, excluded samples are NA.\n")
    (run_dir / "state" / "summary.done.json").unlink(missing_ok=True)
    atomic_json(marker, {"status": "complete", "version": VERSION, "isolates": len(samples), "loci": len(loci),
        "variants": len(variant_rows), "categories": len(category_rows), "signature": signature,
        "outputs": [str(p.relative_to(out)) for folder in (tables, figures, alignments) for p in sorted(folder.rglob("*")) if p.is_file()]})


def scan_outputs_present(out: Path, excluded: bool, cfg: dict) -> bool:
    path = out / "loci.json"
    if not path.is_file():
        return False
    result = load_json(path)
    if excluded:
        return result["status"] == "excluded"
    if not (out / "amrfinder.tsv").is_file():
        return False
    if any(l.get("evidence") in {"assembly_and_reads", "assembly_read_warning"} for l in result["loci"]):
        if not (out / "base_support.tsv.gz").is_file():
            return False
        if truthy(cfg["RESISTANCE_KEEP_BAM"]) and not (out / "reads.bam").is_file():
            if not ((out / "reads.cram").is_file() and (out / "reads.cram.archive.json").is_file()):
                return False
    return True


def prepare(run_dir: Path) -> None:
    """Refresh tool/database provenance and invalidate only changed isolate inputs."""
    from .workers import context
    cfg, _ = context(run_dir)
    preflight(cfg, root_dir(run_dir) / "provenance")
    tasks = read_tsv(run_dir / "state" / "isolate_tasks.tsv")
    for index, task in enumerate(tasks):
        marker = run_dir / "state" / "resistance_scan" / f"{safe_name(task['isolate_id'])}.done.json"
        if marker.is_file():
            _, excluded, signature = sample_inputs(run_dir, index, cfg)
            out = root_dir(run_dir) / "isolates" / safe_name(task["isolate_id"])
            if load_json(marker).get("signature") != signature or not scan_outputs_present(out, excluded, cfg):
                marker.unlink()


def controller(run_dir: Path, cfg: dict) -> None:
    from .workers import _run_index_stage, _run_single_job
    prepare(run_dir)
    tasks = read_tsv(run_dir / "state" / "isolate_tasks.tsv")
    _run_index_stage(run_dir, cfg, "resistance_scan", list(range(len(tasks))), cfg["RESISTANCE_CPUS"],
        cfg["RESISTANCE_MEM"], cfg["RESISTANCE_TIME"], "CleanGene resistance loci")
    _run_single_job(run_dir, cfg, "resistance_merge", cfg["RESISTANCE_MERGE_CPUS"],
        cfg["RESISTANCE_MERGE_MEM"], cfg["RESISTANCE_MERGE_TIME"], "CleanGene resistance alignments")
