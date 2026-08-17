from __future__ import annotations
import shlex, subprocess
from pathlib import Path

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
