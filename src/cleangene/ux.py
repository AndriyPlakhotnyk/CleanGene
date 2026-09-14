from __future__ import annotations
import itertools, sys, threading
from contextvars import ContextVar
from contextlib import contextmanager
from datetime import datetime
import os

RESET="\033[0m"
BURGUNDY="\033[38;2;128;0;32m"
GREEN="\033[38;2;0;128;96m"
OCHRE="\033[38;2;128;96;0m"
BOLD="\033[1m"
SILVER="\033[38;2;192;192;192m"
WHITE="\033[38;2;255;255;255m"

def color_enabled() -> bool:
    mode=os.environ.get("CLEANGENE_COLOR","auto").strip().lower()
    if mode in {"0","false","no","never","off"}: return False
    return mode in {"1","true","yes","always","on"} or sys.stdout.isatty()

def styled(text: str, *, color: str = "", bold: bool = False) -> str:
    if not color_enabled(): return text
    prefix=(BOLD if bold else "") + color
    return f"{prefix}{text}{RESET}" if prefix else text

def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log_line(text: str) -> str:
    return f"{timestamp()} | {text}"

def welcome(text: str) -> str:
    return styled(text,color=BURGUNDY,bold=True)

def waiting(text: str) -> str:
    return styled(text,color=OCHRE)

def completed(text: str) -> str:
    return styled(text,color=GREEN)

def submitted(text: str) -> str:
    return styled(text,color=GREEN,bold=True)

def clean_gene_banner() -> str:
    inner_width=78
    bur=BURGUNDY if color_enabled() else ""
    reset=RESET if color_enabled() else ""
    logo=(
        "  ____   _                           ____                       ",
        " / ___| | |   ___    __ _   _ __    / ___|   ___   _ __     ___ ",
        "| |     | |  / _ \\  / _` | | '_ \\  | |  _   / _ \\ | '_ \\   / _ \\",
        "| |___  | | |  __/ | (_| | | | | | | |_| | |  __/ | | | | |  __/",
        " \\____| |_|  \\___|  \\__,_| |_| |_|  \\____|  \\___| |_| |_|  \\___|",
    )
    rule="=" * inner_width

    def line(text: str = "", color: str = WHITE, bold: bool = False) -> str:
        left=(inner_width-len(text))//2
        right=inner_width-len(text)-left
        return f"{bur}|{' ' * left}{styled(text,color=color,bold=bold)}{bur}{' ' * right}|{reset}"

    def sword() -> str:
        hilt_l="()xxxxx["
        blade_l="======================>"
        gap="    "
        blade_r="<======================"
        hilt_r="]xxxxx()"
        visible=len(hilt_l)+len(blade_l)+len(gap)+len(blade_r)+len(hilt_r)
        left=(inner_width-visible)//2
        right=inner_width-visible-left
        return (
            f"{bur}|{' ' * left}"
            f"{styled(hilt_l,color=BURGUNDY,bold=True)}"
            f"{styled(blade_l+gap+blade_r,color=SILVER,bold=True)}"
            f"{styled(hilt_r,color=BURGUNDY,bold=True)}"
            f"{bur}{' ' * right}|{reset}"
        )

    rows=[f"{bur}+{rule}+{reset}",line(),sword(),line()]
    rows.extend(line(item,GREEN,True) for item in logo)
    rows.extend([
        line(),
        line("Cleanse thy pangenome, my liege!",BURGUNDY,True),
        line(),
        line("Anno Domini 2026",WHITE,True),
        line(),
        sword(),
        line(),
        f"{bur}+{rule}+{reset}",
    ])
    return "\n".join(rows)

_active_spinner = ContextVar("cleangene_spinner", default=None)


class _SpinnerDisplay:
    def __init__(self, message):
        self.message = message
        self.stream = sys.stdout
        self.tty = self.stream.isatty()
        self.lock = threading.RLock()
        self.stop = threading.Event()
        self.width = 0
        self.worker = None

    def clear(self):
        if self.tty:
            self.stream.write("\r" + " " * self.width + "\r")
            self.width = 0

    def set_message(self, message, announce=True):
        with self.lock:
            self.clear()
            self.message = message
            if not self.tty and announce:
                print(message + " | Started", file=self.stream, flush=True)

    def status(self, message, failed):
        with self.lock:
            self.clear()
            text = message + (" | Failed" if failed else " | Complete")
            print(waiting(text) if failed else completed(text), file=self.stream, flush=True)

    def animate(self):
        for mark in itertools.cycle("|/-\\"):
            with self.lock:
                self.clear()
                text = self.message + " " + mark
                self.width = len(text)
                print(waiting(text), end="", file=self.stream, flush=True)
            if self.stop.wait(0.12): break


@contextmanager
def spinner(message: str):
    """One terminal animation, with nested sample/process labels and completion logs."""
    display = _active_spinner.get()
    owner = display is None
    if owner:
        display = _SpinnerDisplay(message)
        token = _active_spinner.set(display)
    previous = display.message
    display.set_message(message)
    if owner and display.tty:
        display.worker = threading.Thread(target=display.animate, daemon=True)
        display.worker.start()
    failed = False
    try:
        yield
    except BaseException:
        failed = True
        raise
    finally:
        if owner:
            display.stop.set()
            if display.worker: display.worker.join()
        display.status(message, failed)
        if owner: _active_spinner.reset(token)
        else: display.set_message(previous, announce=False)
