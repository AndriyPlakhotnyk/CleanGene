from __future__ import annotations
import itertools, sys, threading
from contextlib import contextmanager

@contextmanager
def spinner(message: str):
    stop=threading.Event()
    if not sys.stdout.isatty():
        print(message,flush=True)
        yield
        return
    def animate() -> None:
        for mark in itertools.cycle("|/-\\"):
            print(f"\r{message} {mark}",end="",flush=True)
            if stop.wait(0.12): break
    worker=threading.Thread(target=animate,daemon=True); worker.start()
    failed=False
    try: yield
    except BaseException:
        failed=True
        raise
    finally:
        stop.set(); worker.join(); print(f"\r{message} {'failed' if failed else 'done'}",flush=True)
