#!/usr/bin/env bash
#
# run2.sh
#
# Alternative to run.sh for analogs that only have 100 ns trajectories.
# For each analog:
#   1) extract_metadata.py  -> per-window z(t) time series + WHAM metadata
#      for the last 50 ns (section51-section100) in <out-dir>/<analog>/
#   2) get_pmf2.py          -> WHAM per leaflet, averaged with the standard
#      error estimated from the difference between the two leaflets
#      (solute 1, +z vs. solute 2, -z), using this single 50 ns block.
#
# See run.sh for the primary pipeline (150 ns trajectories, three
# independent 50 ns blocks, block-to-block standard error) using get_pmf.py.
#
# Usage:
#   bash run2.sh                 # all analogs found in ../runs/
#   bash run2.sh scc              # one analog
#   bash run2.sh scc scd scf      # several
#
# Optional environment overrides:
#   PYTHON       python interpreter   (default: python3)
#   WHAM_BIN     Grossfield wham bin  (default: wham)
#   EQUIL_PS     ps to discard/window (default: 0; --min-section already
#                                       starts past equilibration)
#   NBOOT        WHAM bootstrap rounds (default: 0)
#   MIN_SECTION  lowest section index to include  (default: 51)
#   MAX_SECTION  highest section index to include (default: 100)

set -euo pipefail

#HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
HERE="./"
ROOT="/media/bories/Backup_large/bories/Documents/Travail/umbrella/"

RUNS="${RUNS:-$ROOT/runs}"
OUT_ROOT="${OUT_ROOT:-$HERE/pmf2}"
PYTHON="${PYTHON:-python3}"
WHAM_BIN="${WHAM_BIN:-wham}"
EQUIL_PS="${EQUIL_PS:-0}"
NBOOT="${NBOOT:-0}"
MIN_SECTION="${MIN_SECTION:-51}"
MAX_SECTION="${MAX_SECTION:-100}"

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

    "$PYTHON" "$HERE/extract_metadata.py"     \
        --analog       "$a"                   \
        --runs-dir     "$RUNS"                \
        --out-dir      "$OUT_DIR"             \
        --equil-ps     "$EQUIL_PS"            \
        --min-section  "$MIN_SECTION"         \
        --max-section  "$MAX_SECTION"

    "$PYTHON" "$HERE/get_pmf2.py"          \
        --analog   "$a"                   \
        --in-dir   "$OUT_DIR"             \
        --out-dir  "$OUT_DIR"             \
        --wham-bin "$WHAM_BIN"            \
        --nboot    "$NBOOT"
done

echo
echo "Done. PMFs in $OUT_ROOT/<analog>/pmf-us-<analog>.dat"
