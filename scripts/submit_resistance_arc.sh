#!/usr/bin/env bash
# Run from an ARC login node. CleanGene submits a controller and bounded Slurm
# arrays; this launcher does not perform assembly/analysis on the login node.
set -Eeuo pipefail
if (( $# < 2 )); then
    echo "Usage: bash scripts/submit_resistance_arc.sh MANIFEST.tsv CONFIG.env [extra cleangene run options]" >&2
    echo "Activate the cleangene Conda environment first; set ARC account/partition and AMRFINDER_DB in CONFIG.env." >&2
    exit 2
fi
manifest=$1
config=$2
shift 2
project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$project_root/src${PYTHONPATH:+:$PYTHONPATH}"
exec python -m cleangene run \
    --manifest "$manifest" \
    --config "$config" \
    --analysis-root "$project_root" \
    --profile slurm \
    --assembler shovill \
    --skip-downsampling \
    --ignore-checkm2 \
    --resistance-operon \
    "$@"
