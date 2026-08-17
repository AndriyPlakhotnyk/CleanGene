from __future__ import annotations
import shlex, subprocess
from pathlib import Path
from typing import Callable

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

def array_chunks(indices: list[int] | range, chunk_size: int, max_parallel: str) -> list[str]:
    vals=list(indices)
    if chunk_size < 1: raise ValueError("SLURM_ARRAY_CHUNK_SIZE must be >= 1")
    chunks=[]
    for i in range(0,len(vals),chunk_size):
        part=vals[i:i+chunk_size]
        contiguous=part==list(range(part[0],part[-1]+1))
        spec=f"{part[0]}-{part[-1]}" if contiguous and len(part)>1 else ",".join(map(str,part))
        chunks.append(f"{spec}%{max_parallel}")
    return chunks

def submit_chunked_arrays(build_cmd: Callable[[str, str | None], list[str]], arrays: list[str], dry_run: bool, max_outstanding: int = 1, initial_dependency: str | None = None) -> list[str]:
    if max_outstanding < 1: raise ValueError("SLURM_MAX_OUTSTANDING_CHUNKS must be >= 1")
    submitted=[]; wave_dep=initial_dependency
    for i in range(0,len(arrays),max_outstanding):
        wave=[]
        for array in arrays[i:i+max_outstanding]:
            wave.append(submit(build_cmd(array,wave_dep),dry_run))
        submitted.extend(wave)
        wave_dep=":".join(wave)
    return submitted
