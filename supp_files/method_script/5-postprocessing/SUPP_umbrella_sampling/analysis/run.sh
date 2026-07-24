#!/usr/bin/env bash
#
# run.sh
#
# Orchestrate PMF extraction for every analog under ../runs/.
# For each analog:
#   1) extract_metadata.py  -> per-window z(t) time series + WHAM metadata
#   2) get_pmf.py           -> Grossfield WHAM x2 + Allen-2007 4-sample PMF
#
# Usage:
#   bash run.sh                 # all analogs found in ../runs/
#   bash run.sh scc             # one analog
#   bash run.sh scc scd scf     # several
#
# Optional environment overrides:
#   PYTHON       python interpreter   (default: python3)
#   WHAM_BIN     Grossfield wham bin  (default: wham)
#   EQUIL_PS     ps to discard/window (default: 5000)
#   NBOOT        WHAM bootstrap rounds (default: 0)
#   MIN_SECTION  lowest section index to include  (default: 1)
#   MAX_SECTION  highest section index to include (default: no limit)

set -euo pipefail

#HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
HERE="./"
ROOT="../"

RUNS="${RUNS:-$ROOT/runs}"
OUT_ROOT="${OUT_ROOT:-$HERE/pmf}"
PYTHON="${PYTHON:-python3}"
WHAM_BIN="${WHAM_BIN:-wham}"
EQUIL_PS="${EQUIL_PS:-5000}"
NBOOT="${NBOOT:-0}"
MIN_SECTION=100
MAX_SECTION=250

ANALOGS=()
if [ $# -ge 1 ]; then
    ANALOGS=("$@")
else
    for d in "$RUNS"/*/; do
        [ -d "$d" ] || continue
        # only consider directories that actually contain umbrella windows
        if compgen -G "$d/win[0-9][0-9]" > /dev/null; then
            ANALOGS+=("$(basename "$d")")
        fi
    done
fi

if [ ${#ANALOGS[@]} -eq 0 ]; then
    echo "No analogs found under $RUNS"
    exit 0
fi

mkdir -p "$OUT_ROOT"

for a in "${ANALOGS[@]}"; do
    echo "================ $a ================"
    OUT_DIR="$OUT_ROOT/$a"
    mkdir -p "$OUT_DIR"

    EXTRA_EXTRACT_ARGS=(--min-section "$MIN_SECTION")
    if [ -n "$MAX_SECTION" ]; then
        EXTRA_EXTRACT_ARGS+=(--max-section "$MAX_SECTION")
    fi

    "$PYTHON" "$HERE/extract_metadata.py" \
        --analog   "$a"                   \
        --runs-dir "$RUNS"                \
        --out-dir  "$OUT_DIR"             \
        --equil-ps "$EQUIL_PS"            \
        "${EXTRA_EXTRACT_ARGS[@]}"

    "$PYTHON" "$HERE/get_pmf.py"          \
        --analog   "$a"                   \
        --in-dir   "$OUT_DIR"             \
        --out-dir  "$OUT_DIR"             \
        --wham-bin "$WHAM_BIN"            \
        --nboot    "$NBOOT"
done

echo
echo "Done. PMFs in $OUT_ROOT/<analog>/pmf-us-<analog>.dat"
