#!/usr/bin/env bash
#
# run.sh
#
# Orchestrate PMF extraction for every analog under ../runs/.
# For each analog, and for each of three non-overlapping 50 ns blocks
# (section101-150, section151-200, section201-250):
#   1) extract_metadata.py  -> per-window z(t) time series + WHAM metadata
#      for that block, in <out-dir>/<analog>/<block>/
# then, once per analog:
#   2) get_pmf.py           -> WHAM per block (leaflets folded together),
#      averaged over the three blocks. The reported standard error comes
#      from the block-to-block spread, not from the two leaflets.
#
# See run2.sh for the alternative single-block (100 ns trajectory, last
# 50 ns) pipeline using get_pmf2.py, where the SE instead comes from the
# difference between the two leaflets.
#
# Usage:
#   bash run.sh                 # all analogs found in ../runs/
#   bash run.sh scc              # one analog
#   bash run.sh scc scd scf      # several
#
# Optional environment overrides:
#   PYTHON       python interpreter   (default: python3)
#   WHAM_BACKEND python or grossfield (default: python)
#   WHAM_BIN     Grossfield wham bin  (default: wham; used only when
#                                      WHAM_BACKEND=grossfield)
#   EQUIL_PS     ps to discard/block  (default: 0; blocks already start past
#                                       equilibration)
#   NBOOT        WHAM bootstrap rounds (default: 0)

set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT="${ROOT:-/media/bories/Backup_large/bories/Documents/Travail/umbrella}"

RUNS="${RUNS:-$ROOT/runs}"
OUT_ROOT="${OUT_ROOT:-$HERE/pmf}"
PYTHON="${PYTHON:-python3}"
WHAM_BIN="${WHAM_BIN:-wham}"
WHAM_BACKEND="${WHAM_BACKEND:-python}"
EQUIL_PS="${EQUIL_PS:-0}"
NBOOT="${NBOOT:-0}"

# Three contiguous, non-overlapping 50 ns blocks: name -> "min-section max-section"
BLOCK_NAMES=(block1 block2 block3)
BLOCK_RANGES=("101 150" "151 200" "201 250")

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

    for i in "${!BLOCK_NAMES[@]}"; do
        block="${BLOCK_NAMES[$i]}"
        read -r min_section max_section <<< "${BLOCK_RANGES[$i]}"
        echo "-- $a / $block (section$min_section-section$max_section) --"

        "$PYTHON" "$HERE/extract_metadata.py"     \
            --analog       "$a"                   \
            --runs-dir     "$RUNS"                \
            --out-dir      "$OUT_DIR/$block"       \
            --equil-ps     "$EQUIL_PS"             \
            --min-section  "$min_section"          \
            --max-section  "$max_section"
    done

    "$PYTHON" "$HERE/get_pmf.py"          \
        --analog   "$a"                   \
        --in-dir   "$OUT_DIR"             \
        --out-dir  "$OUT_DIR"             \
        --blocks   "${BLOCK_NAMES[@]}"     \
        --backend  "$WHAM_BACKEND"         \
        --wham-bin "$WHAM_BIN"            \
        --nboot    "$NBOOT"
done

echo
echo "Done. PMFs in $OUT_ROOT/<analog>/pmf-us-<analog>.dat"

