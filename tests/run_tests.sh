#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PYTHONPATH="$ROOT/src" python3 -m unittest discover -s "$ROOT/tests" -v
