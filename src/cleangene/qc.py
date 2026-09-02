from __future__ import annotations

import csv, gzip, json, math, shutil
from pathlib import Path
from typing import Iterable

from .util import read_tsv, write_tsv

THRESHOLD_DEFAULTS = {
    "qc_max_contigs_pass": "300",
    "qc_max_contigs_fail": "1000",
    "qc_min_n50_pass": "25000",
    "qc_min_n50_fail": "5000",
    "qc_min_coverage_pass": "20",
    "qc_min_coverage_fail": "10",
    "qc_min_read_length_pass": "120",
    "qc_min_read_length_fail": "",
    "qc_min_mean_base_quality_pass": "30",
    "qc_min_mean_base_quality_fail": "",
    "qc_min_completeness_pass": "90",
    "qc_min_completeness_fail": "80",
    "qc_max_checkm2_contamination_pass": "5",
    "qc_max_checkm2_contamination_fail": "10",
    "qc_max_kraken_contamination_fail": "5",
}
THRESHOLD_COLUMNS = tuple(THRESHOLD_DEFAULTS)
ENV_KEYS = {key:key.upper() for key in THRESHOLD_COLUMNS}
MINIMUM_PAIRS = (
    ("qc_min_n50_pass","qc_min_n50_fail"),
    ("qc_min_coverage_pass","qc_min_coverage_fail"),
    ("qc_min_read_length_pass","qc_min_read_length_fail"),
    ("qc_min_mean_base_quality_pass","qc_min_mean_base_quality_fail"),
    ("qc_min_completeness_pass","qc_min_completeness_fail"),
)
MAXIMUM_PAIRS = (
    ("qc_max_contigs_pass","qc_max_contigs_fail"),
    ("qc_max_checkm2_contamination_pass","qc_max_checkm2_contamination_fail"),
)
QC_OUTPUT_FIELDS = (
    "PASS/FAIL","Notes","trimmed_read_length","mean_base_quality",
    "sequencing_coverage","checkm2_completeness","checkm2_contamination",
    "qc_profile_source",
)

def normalize_organism(value: str | None) -> str:
    return " ".join((value or "").casefold().split())

def _number(value: object, label: str, *, blank: bool = False) -> float | None:
    text=str(value if value is not None else "").strip()
    if not text:
        if blank: return None
        raise SystemExit(f"QC threshold {label} must be numeric")
    try: result=float(text)
    except ValueError: raise SystemExit(f"QC threshold {label} must be numeric: {text}")
    if not math.isfinite(result): raise SystemExit(f"QC threshold {label} must be finite: {text}")
    if result < 0: raise SystemExit(f"QC threshold {label} cannot be negative: {text}")
    return result

def validate_thresholds(values: dict[str,object], label: str = "resolved") -> dict[str,float | None]:
    parsed={key:_number(values.get(key,""),f"{label}.{key}",blank=key in {"qc_min_read_length_fail","qc_min_mean_base_quality_fail"}) for key in THRESHOLD_COLUMNS}
    for pass_key,fail_key in MINIMUM_PAIRS:
        fail=parsed[fail_key]
        if fail is not None and parsed[pass_key] < fail:
            raise SystemExit(f"Invalid QC PASS/FAIL ordering for {label}: {pass_key} must be >= {fail_key}")
    for pass_key,fail_key in MAXIMUM_PAIRS:
        if parsed[pass_key] > parsed[fail_key]:
            raise SystemExit(f"Invalid QC PASS/FAIL ordering for {label}: {pass_key} must be <= {fail_key}")
    return parsed

def _profile_rows(path: Path | None) -> list[dict[str,str]]:
    if not path: return []
    if not path.is_file(): raise SystemExit(f"QC profile file not found: {path}")
    rows=read_tsv(path)
    if not rows:
        header=path.read_text(errors="replace").splitlines()[0].split("\t") if path.read_text(errors="replace").splitlines() else []
    else: header=list(rows[0])
    required={"scope_type","scope_value"}
    if not required.issubset(header): raise SystemExit("QC profile TSV requires scope_type and scope_value columns")
    unknown=set(header)-required-set(THRESHOLD_COLUMNS)
    if unknown: raise SystemExit("Unknown QC profile columns: " + ", ".join(sorted(unknown)))
    seen=set()
    for number,row in enumerate(rows,2):
        scope=row.get("scope_type","").strip(); value=row.get("scope_value","").strip()
        if scope not in {"organism","group_id"}: raise SystemExit(f"Invalid QC profile scope_type on row {number}: {scope}")
        if not value: raise SystemExit(f"QC profile scope_value is blank on row {number}")
        identity=(scope,normalize_organism(value) if scope=="organism" else value)
        if identity in seen: raise SystemExit(f"Duplicate QC profile for {scope}={value}")
        seen.add(identity)
        for key in THRESHOLD_COLUMNS:
            if row.get(key,"").strip(): _number(row[key],f"profile row {number}.{key}")
    return rows

def resolve_threshold_rows(rows: list[dict[str,str]], cfg: dict[str,str], profile_path: Path | None = None) -> list[dict[str,object]]:
    profiles=_profile_rows(profile_path)
    profiles_by_organism={normalize_organism(p["scope_value"]):p for p in profiles if p["scope_type"]=="organism"}
    profiles_by_group={p["scope_value"].strip():p for p in profiles if p["scope_type"]=="group_id"}
    global_values={key:str(cfg.get(ENV_KEYS[key],THRESHOLD_DEFAULTS[key])).strip() for key in THRESHOLD_COLUMNS}
    validate_thresholds(global_values,"global")
    result=[]
    for row in rows:
        values=dict(global_values); sources=[]
        organism=normalize_organism(row.get("organism","")); group=row.get("group_id","").strip()
        for match,source in ((profiles_by_organism.get(organism) if organism else None,f"organism:{row.get('organism','')}"),
                             (profiles_by_group.get(group) if group else None,f"group_id:{group}")):
            if match:
                changed=False
                for key in THRESHOLD_COLUMNS:
                    if match.get(key,"").strip(): values[key]=match[key].strip(); changed=True
                if changed: sources.append(source)
        manifest_changed=False
        for key in THRESHOLD_COLUMNS:
            if row.get(key,"").strip():
                _number(row[key],f"manifest {row['isolate_id']}.{key}")
                values[key]=row[key].strip(); manifest_changed=True
        if manifest_changed: sources.append(f"manifest:{row['isolate_id']}")
        validate_thresholds(values,row["isolate_id"])
        source=sources[-1] if sources else "global"
        result.append({"isolate_id":row["isolate_id"],"qc_profile_source":source,**values})
    return result

def prepare_qc_provenance(run_dir: Path, rows: list[dict[str,str]], cfg: dict[str,str]) -> dict[str,str]:
    cfg=dict(cfg)
    if cfg.get("CHECKM2_DB","").strip(): cfg["CHECKM2_DB"]=str(Path(cfg["CHECKM2_DB"]).expanduser().resolve())
    configured=cfg.get("QC_PROFILE_FILE","").strip(); profile=Path(configured).expanduser().resolve() if configured else None
    if profile:
        if not profile.is_file(): raise SystemExit(f"QC profile file not found: {profile}")
        copied=run_dir/"provenance"/"qc_profile.tsv"; copied.parent.mkdir(parents=True,exist_ok=True)
        if profile != copied.resolve(): shutil.copy2(profile,copied)
        profile=copied; cfg["QC_PROFILE_FILE"]=str(copied)
    resolved=resolve_threshold_rows(rows,cfg,profile)
    write_tsv(run_dir/"provenance"/"qc_thresholds.tsv",["isolate_id","qc_profile_source",*THRESHOLD_COLUMNS],resolved)
    return cfg

def ensure_qc_provenance(run_dir: Path, rows: list[dict[str,str]], cfg: dict[str,str]) -> dict[str,str]:
    target=run_dir/"provenance"/"qc_thresholds.tsv"
    if target.is_file(): return cfg
    return prepare_qc_provenance(run_dir,rows,cfg)

def isolate_thresholds(run_dir: Path, isolate_id: str, cfg: dict[str,str], row: dict[str,str]) -> tuple[dict[str,float | None],str]:
    path=run_dir/"provenance"/"qc_thresholds.tsv"
    matches=[item for item in read_tsv(path) if item.get("isolate_id")==isolate_id] if path.is_file() else resolve_threshold_rows([row],cfg,None)
    if len(matches)!=1: raise SystemExit(f"Resolved QC thresholds missing or duplicated for isolate {isolate_id}")
    return validate_thresholds(matches[0],isolate_id),str(matches[0].get("qc_profile_source","global"))

def _open_fastq(path: Path):
    with path.open("rb") as handle: gz=handle.read(2)==b"\x1f\x8b"
    return gzip.open(path,"rt",errors="replace") if gz else path.open(errors="replace")

def _scan_fastq(path: Path) -> tuple[int,int,int]:
    reads=bases=quality_sum=0
    with _open_fastq(path) as handle:
        while True:
            header=handle.readline()
            if not header: break
            sequence=handle.readline().rstrip("\r\n"); plus=handle.readline(); quality=handle.readline().rstrip("\r\n")
            if not sequence or not plus or len(sequence)!=len(quality): raise SystemExit(f"Malformed FASTQ while computing QC metrics: {path}")
            reads += 1; bases += len(sequence); quality_sum += sum(ord(char)-33 for char in quality)
    if reads==0 or bases==0: raise SystemExit(f"FASTQ contains no reads: {path}")
    return reads,bases,quality_sum

def read_metrics(r1: Path, r2: Path, fastp_json: Path | None = None) -> dict[str,float]:
    r1_reads,r1_bases,r1_quality=_scan_fastq(r1); r2_reads,r2_bases,r2_quality=_scan_fastq(r2)
    total_bases=r1_bases+r2_bases; mean1=r1_bases/r1_reads; mean2=r2_bases/r2_reads
    if fastp_json and fastp_json.is_file():
        try:
            after=json.loads(fastp_json.read_text()).get("summary",{}).get("after_filtering",{})
            candidate_total=int(after["total_bases"]); candidate1=float(after["read1_mean_length"]); candidate2=float(after["read2_mean_length"])
            if candidate_total>0 and candidate1>0 and candidate2>0:
                total_bases=candidate_total; mean1=candidate1; mean2=candidate2
        except (OSError,ValueError,TypeError,KeyError,json.JSONDecodeError): pass
    return {"trimmed_read_length":min(mean1,mean2),"mean_base_quality":(r1_quality+r2_quality)/(r1_bases+r2_bases),"total_bases":float(total_bases)}

def parse_checkm2_report(path: Path) -> tuple[float,float]:
    if not path.is_file(): raise ValueError(f"CheckM2 report not found: {path}")
    with path.open(newline="",errors="replace") as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    if not rows: raise ValueError(f"CheckM2 report contains no genome rows: {path}")
    normalized={key.casefold().strip():value for key,value in rows[0].items() if key}
    try: return float(normalized["completeness"]),float(normalized["contamination"])
    except (KeyError,TypeError,ValueError) as error: raise ValueError(f"Invalid CheckM2 completeness/contamination in {path}") from error

def _display(value: float) -> str:
    number=float(value)
    return str(int(number)) if number.is_integer() else str(number)

def classify_isolate_qc(*, thresholds: dict[str,float | None], expected_organism: str, top_species: str,
                        kraken_contamination: float | None, read_length: float | None,
                        mean_quality: float | None, coverage: float | None, contigs: float | None,
                        n50: float | None, completeness: float | None, checkm2_contamination: float | None,
                        checkm2_mode: str, internal_pangenome: bool, external_pangenome: bool,
                        assembly_present: bool, gff_present: bool | None, user_exclusion: bool = False,
                        warnings: Iterable[tuple[str,str]] = (), errors: Iterable[tuple[str,str]] = ()) -> dict[str,str]:
    notes=[]; info_notes=[]; codes=[]
    def add(level: str,code: str,message: str) -> None:
        notes.append((level,message))
        if level=="FAIL": codes.append(code)
    def info(message: str) -> None:
        info_notes.append(message)
    if user_exclusion:
        add("FAIL","user_excluded","isolate was explicitly excluded by the user")
        return {"PASS/FAIL":"FAIL","Notes":"FAIL: isolate was explicitly excluded by the user","excluded":"1","reason":"user_excluded"}
    expected=normalize_organism(expected_organism); observed=normalize_organism(top_species)
    if expected:
        if not observed: add("WARNING","taxonomy_unavailable",f"expected organism '{expected_organism}' was supplied but Kraken classification was unavailable")
        elif observed!=expected: add("FAIL","organism_mismatch",f"Kraken top species '{top_species}' does not match expected organism '{expected_organism}'")
    if kraken_contamination is not None and kraken_contamination>thresholds["qc_max_kraken_contamination_fail"]:
        add("FAIL","kraken_contamination_high",f"Kraken contamination={_display(kraken_contamination)}% exceeds fail maximum {_display(thresholds['qc_max_kraken_contamination_fail'])}%")
    def minimum(value: float | None, pass_key: str, fail_key: str, label: str, unit: str, unavailable: str | None = None) -> None:
        code_label=label.casefold()
        if value is None:
            if unavailable: add("WARNING",f"{code_label}_unavailable",unavailable)
            return
        fail=thresholds[fail_key]
        if fail is not None and value<fail: add("FAIL",f"{code_label}_low",f"{label}={_display(value)}{unit} is below fail minimum {_display(fail)}{unit}")
        elif value<thresholds[pass_key]: add("WARNING",f"{code_label}_warning",f"{label}={_display(value)}{unit} is below pass minimum {_display(thresholds[pass_key])}{unit}")
    def maximum(value: float | None, pass_key: str, fail_key: str, label: str, unit: str = "") -> None:
        code_label=label.casefold()
        if value is None: return
        if value>thresholds[fail_key]: add("FAIL",f"{code_label}_high",f"{label}={_display(value)}{unit} exceeds fail maximum {_display(thresholds[fail_key])}{unit}")
        elif value>thresholds[pass_key]: add("WARNING",f"{code_label}_warning",f"{label}={_display(value)}{unit} exceeds pass maximum {_display(thresholds[pass_key])}{unit}")
    minimum(read_length,"qc_min_read_length_pass","qc_min_read_length_fail","read_length"," bp","post-processing read length was unavailable")
    minimum(mean_quality,"qc_min_mean_base_quality_pass","qc_min_mean_base_quality_fail","mean_base_quality","","post-processing mean base quality was unavailable")
    minimum(coverage,"qc_min_coverage_pass","qc_min_coverage_fail","coverage","x","sequencing coverage was unavailable")
    maximum(contigs,"qc_max_contigs_pass","qc_max_contigs_fail","contigs")
    minimum(n50,"qc_min_n50_pass","qc_min_n50_fail","N50"," bp")
    if checkm2_mode=="off": info("CheckM2 was not evaluated because it was explicitly disabled")
    else:
        minimum(completeness,"qc_min_completeness_pass","qc_min_completeness_fail","completeness","%","CheckM2 completeness was unavailable")
        maximum(checkm2_contamination,"qc_max_checkm2_contamination_pass","qc_max_checkm2_contamination_fail","checkm2_contamination","%")
    for code,message in warnings: add("WARNING",code,message)
    for code,message in errors: add("FAIL",code,message)
    if external_pangenome and not assembly_present:
        add("WARNING","external_pangenome_assembly_unavailable","assembly metrics and Prokka GFF were not evaluated because an external pangenome was supplied without an assembly")
    elif external_pangenome:
        add("WARNING","external_pangenome_gff_not_evaluated","Prokka GFF was not evaluated because an external pangenome was supplied")
    elif internal_pangenome:
        if not assembly_present: add("FAIL","assembly_missing","assembly was not produced for the internal pangenome")
        elif gff_present is False: add("FAIL","gff_missing","Prokka GFF was not produced for the internal pangenome")
    level="FAIL" if any(item[0]=="FAIL" for item in notes) else "WARNING" if notes else "PASS"
    rendered_notes=[f"{severity}: {message}" for severity,message in notes]
    rendered_notes.extend(f"INFO: {message}" for message in info_notes)
    rendered="; ".join(rendered_notes) if rendered_notes else "All evaluated QC criteria passed"
    return {"PASS/FAIL":level,"Notes":rendered,"excluded":"1" if level=="FAIL" else "0","reason":";".join(codes)}

def qc_value(value: float | int | None) -> str:
    if value is None: return ""
    return _display(float(value))
