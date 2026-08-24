from __future__ import annotations
import itertools, sys, threading
from contextlib import contextmanager

RESET="\033[0m"
BURGUNDY="\033[38;2;128;0;32m"
GREEN="\033[38;2;0;128;96m"
OCHRE="\033[38;2;128;96;0m"
BOLD="\033[1m"

def styled(text: str, *, color: str = "", bold: bool = False) -> str:
    prefix=(BOLD if bold else "") + color
    return f"{prefix}{text}{RESET}" if prefix else text

def welcome(text: str) -> str:
    return styled(text,color=BURGUNDY,bold=True)

def waiting(text: str) -> str:
    return styled(text,color=OCHRE)

def completed(text: str) -> str:
    return styled(text,color=GREEN)

def submitted(text: str) -> str:
    return styled(text,color=GREEN,bold=True)

@contextmanager
def spinner(message: str):
    stop=threading.Event()
    display=waiting(message)
    if not sys.stdout.isatty():
        print(display,flush=True)
        yield
        return
    def animate() -> None:
        for mark in itertools.cycle("|/-\\"):
            print(f"\r{display} {waiting(mark)}",end="",flush=True)
            if stop.wait(0.12): break
    worker=threading.Thread(target=animate,daemon=True); worker.start()
    failed=False
    try: yield
    except BaseException:
        failed=True
        raise
    finally:
        stop.set(); worker.join(); status="failed" if failed else "done"; print(f"\r{display} {waiting(status)}",flush=True)
