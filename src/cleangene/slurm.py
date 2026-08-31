from __future__ import annotations
import shlex, subprocess
import os, time
from collections import Counter
from pathlib import Path
from .ux import log_line, waiting

def sbatch_cmd(*, name: str, wrap: str, cpus: str, mem: str, time: str, array: str | None = None, dependency: str | None = None, account: str = "", partition: str = "", log: Path | None = None) -> list[str]:
    cmd=["sbatch","--parsable","--job-name",name,"--cpus-per-task",cpus,"--mem",mem,"--time",time]
    if array: cmd += ["--array",array]
    if dependency: cmd += ["--dependency",f"afterok:{dependency}","--kill-on-invalid-dep=yes"]
    if account: cmd += ["--account",account]
    if partition: cmd += ["--partition",partition]
    if log: cmd += ["--output",str(log),"--error",str(log)]
    cmd += ["--wrap",wrap]
    return cmd

def submit(cmd: list[str], dry_run: bool) -> str:
    if dry_run:
        print(" ".join(shlex.quote(x) for x in cmd)); return "DRYRUN"
    result=subprocess.run(cmd,capture_output=True,text=True)
    if result.returncode:
        rendered=" ".join(shlex.quote(x) for x in cmd)
        stderr=result.stderr.strip() or "<no stderr>"
        stdout=result.stdout.strip()
        message=f"sbatch failed with exit status {result.returncode}\ncommand: {rendered}\nstderr: {stderr}"
        if stdout:
            message += f"\nstdout: {stdout}"
        raise RuntimeError(message)
    return result.stdout.strip().split(";")[0]

def array_task_count(array: str | None) -> int:
    if not array: return 1
    spec=array.split("%",1)[0]; total=0
    for part in spec.split(","):
        if "-" in part:
            a,b=part.split("-",1); total += int(b)-int(a)+1
        elif part:
            total += 1
    return total

def user_job_count(user: str | None = None) -> int:
    return user_queue_snapshot(user)["total"]

def user_queue_snapshot(user: str | None = None) -> dict[str,object]:
    user=user or os.environ.get("USER","")
    result=subprocess.run(["squeue","-r","-h","-u",user,"-t","R,PD","-o","%F|%K|%T|%j|%C|%o"],capture_output=True,text=True,check=True)
    jobs: dict[str,Counter[str]]= {}
    entries=[]
    total=0
    for line in result.stdout.splitlines():
        if not line.strip(): continue
        fields=line.strip().split("|",5)
        job_id=fields[0].strip()
        task_id=fields[1].strip() if len(fields)>1 else ""
        state=fields[2].strip().upper() if len(fields)>2 else "UNKNOWN"
        state={"R":"RUNNING","PD":"PENDING"}.get(state,state or "UNKNOWN")
        name=fields[3].strip() if len(fields)>3 else ""
        cpus=1
        command=""
        if len(fields)>5:
            try: cpus=max(0,int(fields[4].strip() or "0"))
            except ValueError: cpus=0
            command=fields[5].strip()
        elif len(fields)>4:
            command=fields[4].strip()
        jobs.setdefault(job_id,Counter())[state] += 1
        entries.append({"job_id":job_id,"task_id":task_id,"state":state,"name":name,"cpus":cpus,"command":command})
        total += 1
    return {"total":total,"jobs":jobs,"entries":entries}

def active_cleangene_jobs_for_run(run_dir: Path, user: str | None = None) -> list[dict[str,object]]:
    run_arg=str(run_dir.expanduser().resolve())
    try:
        snapshot=user_queue_snapshot(user)
    except (OSError, subprocess.SubprocessError):
        return []
    matches=[]
    for entry in snapshot.get("entries",[]):
        command=str(entry.get("command",""))
        name=str(entry.get("name",""))
        if run_arg in command and name.startswith("cg-"):
            matches.append(entry)
    return matches

def cancel_jobs(job_ids: list[str]) -> None:
    unique=sorted({jid for jid in job_ids if jid and jid!="DRYRUN"})
    if unique:
        subprocess.run(["scancel", *unique], check=True)

def available_slots(limit: int, headroom: int, current: int) -> int:
    return max(0, limit - headroom - current)

def wait_for_capacity(needed: int, cfg: dict[str,str], *, label: str = "") -> int:
    limit=int(cfg["SLURM_USER_JOB_LIMIT"]); headroom=int(cfg["SLURM_JOB_HEADROOM"]); poll=int(cfg["SLURM_POLL_SECONDS"])
    while True:
        current=user_job_count()
        avail=available_slots(limit,headroom,current)
        print(waiting(log_line(f"step={label} | user_jobs={current}/{limit} | available_slots={avail} | waiting_for_capacity | needed={needed}")), flush=True)
        if avail >= needed: return avail
        time.sleep(poll)

def submit_with_qos_retry(cmd: list[str], cfg: dict[str,str], task_count: int, label: str = "") -> str:
    poll=int(cfg["SLURM_POLL_SECONDS"])
    while True:
        wait_for_capacity(task_count,cfg,label=label)
        try:
            return submit(cmd,False)
        except RuntimeError as e:
            if "QOSMaxSubmitJobPerUserLimit" not in str(e) and "violates accounting/QOS policy" not in str(e):
                raise
            print(str(e), flush=True)
            time.sleep(poll)

def job_active(job_ids: list[str]) -> set[str]:
    ids=[j for j in job_ids if j and j!="DRYRUN"]
    if not ids: return set()
    result=subprocess.run(["squeue","-h","-j",",".join(ids),"-o","%A"],capture_output=True,text=True)
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}

def assert_jobs_succeeded(job_ids: list[str], details: str = "") -> None:
    ids=[j for j in job_ids if j and j!="DRYRUN"]
    if not ids: return
    result=subprocess.run(["sacct","-n","-X","-j",",".join(ids),"--format=JobIDRaw,State","-P"],capture_output=True,text=True)
    if result.returncode: return
    bad=[s.strip() for s in result.stdout.splitlines() if s.strip() and not s.split("|")[-1].startswith(("COMPLETED","COMPLETING"))]
    if bad: raise RuntimeError("SLURM job failure detected: " + ", ".join(bad[:5]) + (f" | {details}" if details else ""))
