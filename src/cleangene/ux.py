from __future__ import annotations
import itertools, sys, threading
from contextlib import contextmanager

RESET="\033[0m"
BURGUNDY="\033[38;2;128;0;32m"
GREEN="\033[38;2;0;128;96m"
OCHRE="\033[38;2;128;96;0m"
BOLD="\033[1m"
SILVER="\033[38;2;192;192;192m"
WHITE="\033[38;2;255;255;255m"

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

def clean_gene_banner() -> str:
    inner_width=78
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
        return f"{BURGUNDY}|{' ' * left}{styled(text,color=color,bold=bold)}{BURGUNDY}{' ' * right}|{RESET}"

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
            f"{BURGUNDY}|{' ' * left}"
            f"{styled(hilt_l,color=BURGUNDY,bold=True)}"
            f"{styled(blade_l+gap+blade_r,color=SILVER,bold=True)}"
            f"{styled(hilt_r,color=BURGUNDY,bold=True)}"
            f"{BURGUNDY}{' ' * right}|{RESET}"
        )

    rows=[f"{BURGUNDY}+{rule}+{RESET}",line(),sword(),line()]
    rows.extend(line(item,GREEN,True) for item in logo)
    rows.extend([
        line(),
        line("Cleanse thy pangenome, my liege!",BURGUNDY,True),
        line(),
        line("Anno Domini 2026",WHITE,True),
        line(),
        sword(),
        line(),
        f"{BURGUNDY}+{rule}+{RESET}",
    ])
    return "\n".join(rows)

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
