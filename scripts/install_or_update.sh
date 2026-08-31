#!/usr/bin/env bash
set -euo pipefail

trap 'echo "ERROR: CleanGene installation/update failed." >&2' ERR

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

usage() {
    echo "Usage: bash scripts/install_or_update.sh [--recreate]" >&2
}

RECREATE=false
if [[ ${1:-} == "--recreate" ]]; then
    RECREATE=true
    shift
fi
if [[ $# -ne 0 ]]; then
    usage
    exit 2
fi
if ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: mamba is required. Install Miniforge/Mambaforge, then rerun this script." >&2
    exit 1
fi

environment_exists() {
    local name=$1
    mamba env list | awk -v target="$name" '$1 == target { found=1 } END { exit !found }'
}

install_environment() {
    local name=$1
    local file=$2
    if [[ $RECREATE == true ]]; then
        if environment_exists "$name"; then
            echo "Recreating environment: $name"
            mamba env remove -n "$name" --yes
        else
            echo "Creating environment: $name"
        fi
        mamba env create -f "$file"
    elif environment_exists "$name"; then
        echo "Updating environment: $name"
        mamba env update -n "$name" -f "$file" --prune
    else
        echo "Creating environment: $name"
        mamba env create -f "$file"
    fi
}

install_environment cleangene environment.yml
install_environment cleangene-checkm2 environment.checkm2.yml

echo "Installing this CleanGene checkout"
mamba run -n cleangene python -m pip install -e .

echo "Verifying CleanGene editable checkout"
mamba run -n cleangene python -c '
import cleangene, pathlib, subprocess, sys
root = pathlib.Path.cwd().resolve()
package = pathlib.Path(cleangene.__file__).resolve().parent
expected = root / "src" / "cleangene"
commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
print(f"CleanGene package path: {package}")
print(f"CleanGene expected path: {expected}")
print(f"CleanGene python: {pathlib.Path(sys.executable).resolve()}")
print(f"CleanGene git HEAD: {commit}")
if package != expected:
    raise SystemExit(f"ERROR: cleangene imports from {package}, expected {expected}")
'

echo "Verifying CheckM2 predict CLI"
mamba run -n cleangene python -c '
from cleangene.checkm2 import checkm2_predict_capabilities
from cleangene.tools import executable_version, resolve_checkm2_executable
exe = resolve_checkm2_executable("")
version = executable_version(exe, "CheckM2")
cap = checkm2_predict_capabilities(exe)
print(f"CheckM2 executable: {exe}")
print(f"CheckM2 version: {version}")
print(f"CheckM2 cleanup option: {cap.cleanup_option}")
'

echo "Verifying primary tools"
mamba run -n cleangene sh -c '
    for tool in shovill spades.py prokka panaroo bwa samtools bcftools minimap2 fastp kraken2; do
        command -v "$tool" >/dev/null || { echo "ERROR: missing primary tool: $tool" >&2; exit 1; }
    done
'

if [[ ! -e config/cleangene.arc.local.env ]]; then
    cp config/cleangene.arc.env config/cleangene.arc.local.env
    echo "Created config/cleangene.arc.local.env"
else
    echo "Keeping existing config/cleangene.arc.local.env"
fi

echo "Running CleanGene deployment checks"
mamba run -n cleangene cleangene doctor --config config/cleangene.arc.local.env

echo "CleanGene installation/update completed successfully."
echo "Next command:"
echo "  conda activate cleangene"
