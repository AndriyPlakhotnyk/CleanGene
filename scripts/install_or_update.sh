#!/usr/bin/env bash
set -Eeuo pipefail

trap 'echo "ERROR: CleanGene installation/update failed." >&2' ERR

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

usage() {
    echo "Usage: bash scripts/install_or_update.sh [--recreate] [--env-root DIRECTORY] [--profile local|slurm]" >&2
}

RECREATE=false
ENV_ROOT=""
PROFILE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --recreate) RECREATE=true; shift ;;
        --profile)
            [[ $# -ge 2 && ( "$2" == local || "$2" == slurm ) ]] || { usage; exit 2; }
            PROFILE=$2; shift 2 ;;
        --env-root)
            [[ $# -ge 2 ]] || { usage; exit 2; }
            mkdir -p "$2"
            ENV_ROOT=$(cd "$2" && pwd)
            shift 2 ;;
        *) usage; exit 2 ;;
    esac
done
if ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: mamba is required. Install Miniforge/Mambaforge, then rerun this script." >&2
    exit 1
fi

environment_exists() {
    local name=$1
    if [[ -n "$ENV_ROOT" ]]; then
        [[ -d "$ENV_ROOT/$name/conda-meta" ]]
        return
    fi
    mamba env list | awk -v target="$name" '$1 == target { found=1 } END { exit !found }'
}

install_environment() {
    local name=$1
    local file=$2
    local selector=(-n "$name")
    local create_selector=()
    if [[ -n "$ENV_ROOT" ]]; then
        selector=(-p "$ENV_ROOT/$name")
        create_selector=("${selector[@]}")
    fi
    if [[ $RECREATE == true ]]; then
        if environment_exists "$name"; then
            echo "Recreating environment: $name"
            mamba env remove "${selector[@]}" --yes
        else
            echo "Creating environment: $name"
        fi
        mamba env create "${create_selector[@]}" -f "$file" --yes
    elif environment_exists "$name"; then
        echo "Updating environment: $name"
        mamba env update "${selector[@]}" -f "$file" --prune --yes
    else
        echo "Creating environment: $name"
        mamba env create "${create_selector[@]}" -f "$file" --yes
    fi
}

run_primary() {
    if [[ -n "$ENV_ROOT" ]]; then
        mamba run -p "$ENV_ROOT/cleangene" "$@"
    else
        mamba run -n cleangene "$@"
    fi
}

install_environment cleangene environment.yml
install_environment cleangene-checkm2 environment.checkm2.yml

echo "Installing this CleanGene checkout"
run_primary python -m pip install -e .

echo "Verifying CleanGene editable checkout"
run_primary python -c '
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
run_primary python -c '
from cleangene.checkm2 import bundled_test_genome, checkm2_predict_capabilities
from cleangene.tools import executable_version, resolve_checkm2_executable
exe = resolve_checkm2_executable("")
version = executable_version(exe, "CheckM2")
cap = checkm2_predict_capabilities(exe)
print(f"CheckM2 executable: {exe}")
print(f"CheckM2 version: {version}")
print(f"CheckM2 cleanup option: {cap.cleanup_option}")
print(f"CheckM2 bundled test genome: {bundled_test_genome(exe)}")
'

echo "Verifying primary tools"
run_primary sh -c '
    for tool in shovill spades.py prokka panaroo bwa samtools bcftools minimap2 fastp kraken2 prodigal cd-hit-est cd-hit-est-2d mafft; do
        command -v "$tool" >/dev/null || { echo "ERROR: missing primary tool: $tool" >&2; exit 1; }
    done
'

echo "Verifying Shovill dependencies"
run_primary shovill --check

echo "Verifying Panaroo runtime"
run_primary python -c '
from importlib.metadata import version
from panaroo.__main__ import main
package = "panaroo"
print(f"Panaroo version: {version(package)}")
print("Panaroo Python entry point: READY")
'

if [[ ! -e config/cleangene.arc.local.env ]]; then
    cp config/cleangene.arc.env config/cleangene.arc.local.env
    echo "Created config/cleangene.arc.local.env"
else
    echo "Keeping existing config/cleangene.arc.local.env"
fi

if [[ -z "$PROFILE" ]]; then
    if command -v sbatch >/dev/null 2>&1; then PROFILE=slurm; else PROFILE=local; fi
fi
echo "Running CleanGene deployment checks (profile: $PROFILE)"
run_primary cleangene doctor --config config/cleangene.arc.local.env --profile "$PROFILE"

echo "CleanGene installation/update completed successfully."
echo "Next command:"
if [[ -n "$ENV_ROOT" ]]; then
    echo "  conda activate $ENV_ROOT/cleangene"
else
    echo "  conda activate cleangene"
fi
