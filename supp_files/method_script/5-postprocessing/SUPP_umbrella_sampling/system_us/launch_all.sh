#!/usr/bin/env bash
#
# launch_all.sh
#
# Submit every prepared umbrella window (3 analogs x 38 windows = 114 jobs)
# to the SLURM scheduler. Each job script is self-resubmitting until 30 ns
# of production is reached (see templates/submit.sh).
#
# Usage:
#       bash launch_all.sh                  # submit everything
#       bash launch_all.sh scc              # only one analog
#       bash launch_all.sh scd 05 06 07     # only listed windows of one analog

set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
RUNS_DIR="$HERE/runs"
ANALOGS=(scl)

WINDOWS=()
for i in $(seq 0 37); do WINDOWS+=($(printf "%02d" "$i")); done


if [ $# -ge 1 ]; then
    ANALOGS=("$1"); shift
fi
if [ $# -ge 1 ]; then
    WINDOWS=("$@")
fi

submitted=0
for a in "${ANALOGS[@]}"; do
    for w in "${WINDOWS[@]}"; do
        DIR="$RUNS_DIR/$a/win$w"
        if [ ! -f "$DIR/submit.sh" ]; then
            echo "  ! missing $DIR/submit.sh -- did you run setup_all.sh ?"
            continue
        fi
        if [ -f "$DIR/out/section250.out.gz" ] || [ -f "$DIR/out/section250.out" ]; then
            echo "  - $a/win$w already finished, skipping"
            continue
        fi
        ( cd "$DIR" && sbatch submit.sh )
        submitted=$((submitted + 1))
    done
done

echo
echo "Submitted $submitted jobs."
