#!/usr/bin/env bash
# Update or clone CleanGene on an ARC login node, install the current software,
# and provision a versioned AMRFinderPlus database.
set -Eeuo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  bash scripts/bootstrap_arc_latest.sh REPOSITORY_DIR ENV_ROOT AMRFINDER_DB_DIR [BRANCH]

Example:
  bash scripts/bootstrap_arc_latest.sh \
    /project/def-yourlab/cleangene \
    /project/def-yourlab/conda-envs \
    /project/def-yourlab/databases/amrfinderplus main

Before running, load ARC's Git and Miniforge/Mamba modules if needed. The
repository must be the CleanGene GitHub checkout or an existing clone with the
same origin. ENV_ROOT and AMRFINDER_DB_DIR should be on project/shared storage,
not a compute-node temporary directory.
EOF
    exit 2
}

[[ $# -ge 3 && $# -le 4 ]] || usage
REPO_DIR=$(readlink -f "$1")
ENV_ROOT=$(readlink -f "$2")
AMRFINDER_DB_DIR=$(readlink -f "$3")
BRANCH=${4:-main}
REMOTE_URL=${CLEANGENE_REMOTE:-git@github.com:AndriyPlakhotnyk/CleanGene.git}

command -v git >/dev/null || { echo "ERROR: git is not on PATH; load the ARC Git module." >&2; exit 1; }
command -v mamba >/dev/null || {
    echo "ERROR: mamba is not on PATH; load Miniforge/Mambaforge or set up its shell hook." >&2
    exit 1
}

if [[ -e "$REPO_DIR" && ! -d "$REPO_DIR/.git" ]]; then
    echo "ERROR: $REPO_DIR exists but is not a Git checkout." >&2
    exit 1
fi

if [[ ! -d "$REPO_DIR/.git" ]]; then
    mkdir -p "$(dirname "$REPO_DIR")"
    git clone --origin origin "$REMOTE_URL" "$REPO_DIR"
fi

cd "$REPO_DIR"
git remote get-url origin >/dev/null
if [[ -n "$(git status --porcelain)" ]]; then
    echo "ERROR: checkout has local changes; refusing to overwrite them:" >&2
    git status --short >&2
    echo "Commit/stash those changes, then rerun this script." >&2
    exit 1
fi
git fetch --prune origin
if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
    git checkout "$BRANCH"
else
    git checkout -B "$BRANCH" "origin/$BRANCH"
fi
git pull --ff-only origin "$BRANCH"

mkdir -p "$ENV_ROOT" "$AMRFINDER_DB_DIR"
bash scripts/install_or_update.sh --env-root "$ENV_ROOT" --profile slurm

ENV_PREFIX="$ENV_ROOT/cleangene"
echo "Updating AMRFinderPlus database in $AMRFINDER_DB_DIR"
mamba run -p "$ENV_PREFIX" amrfinder_update --database "$AMRFINDER_DB_DIR"

DB_VERSION=$(find "$AMRFINDER_DB_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V | tail -n 1)
[[ -n "$DB_VERSION" ]] || { echo "ERROR: AMRFinderPlus database update produced no version directory." >&2; exit 1; }
VERSIONED_DB="$AMRFINDER_DB_DIR/$DB_VERSION"

LOCAL_CONFIG="$REPO_DIR/config/cleangene.arc.local.env"
if [[ ! -e "$LOCAL_CONFIG" ]]; then
    cp config/cleangene.arc.env "$LOCAL_CONFIG"
fi
if ! grep -q '^AMRFINDER_DB=' "$LOCAL_CONFIG"; then
    printf '\nAMRFINDER_DB="%s"\n' "$VERSIONED_DB" >> "$LOCAL_CONFIG"
else
    sed -i "s|^AMRFINDER_DB=.*|AMRFINDER_DB=\"$VERSIONED_DB\"|" "$LOCAL_CONFIG"
fi

echo "CleanGene commit: $(git rev-parse --short HEAD)"
echo "AMRFinderPlus database: $VERSIONED_DB"
mamba run -p "$ENV_PREFIX" cleangene doctor --profile slurm --config "$LOCAL_CONFIG"
mamba run -p "$ENV_PREFIX" amrfinder --database "$VERSIONED_DB" --database_version
echo "ARC bootstrap completed. Activate with:"
echo "  conda activate $ENV_PREFIX"
echo "Then edit $LOCAL_CONFIG for SLURM_ACCOUNT, SLURM_PARTITION, and your input paths."
