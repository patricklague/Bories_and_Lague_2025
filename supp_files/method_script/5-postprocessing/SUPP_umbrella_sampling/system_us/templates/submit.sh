#!/bin/bash
#SBATCH --account=def-plague
#SBATCH --job-name=US__ANALOG____WIN__
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=2048
#SBATCH --time=0-24:00          # 24 h per submission; resubmits until NSECTION
#SBATCH --output=slurm-%j.out
#
# Per-window driver for umbrella sampling.
# Mirrors the section-based pattern from membrane/charmm-gui/namd/submit.sh:
#   - 6-step CHARMM-GUI equilibration (step6.1 .. step6.6), step6.1 includes minimize
#   - then N sections of production (NSECTION * 1 ns)
#   - each section produces out/sectionN.{dcd,coor,vel,xsc,xst,colvars.{state,traj}}
#   - script resubmits itself until NSECTION sections exist
#
# Substituted by setup_all.sh: scc, 00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load StdEnv/2023
module load gcc/12.3
module load cuda/12.6
module load namd-multicore/3.0.1

NAMD="namd3 +p$SLURM_CPUS_PER_TASK +idlepoll"
#NAMD="namd3 +p5 +idlepoll"

WIN=__WIN__
ANALOG=__ANALOG__

# Total production = NSECTION * 1 ns  (each section is 500000 steps * 2 fs)
NSECTION=250
SECTIONS_PER_SUBMISSION=50

# --------------------------------------------------------- equilibration (6 steps)
# CHARMM-GUI-style ramp; step6.1 also performs the initial minimization.
for s in 1 2 3 4 5 6; do
    if [ ! -f step6.${s}_eq.restart.vel ]; then
        $NAMD step6.${s}_equilibration.namd > step6.${s}_eq.out
    fi
done

# ------------------------------------------------------------------ production
mkdir -p out

# On first production submission, seed out/restart.* from step6.6
if [ ! -f out/section1.vel.gz ]; then
    cp step6.6_eq.restart.coor out/restart.coor
    cp step6.6_eq.restart.vel  out/restart.vel
    cp step6.6_eq.restart.xsc  out/restart.xsc
    # carry over the colvars bias state so production picks up where step6.6 left off
    if [ -f step6.6_eq.colvars.state ]; then
        cp step6.6_eq.colvars.state out/restart.colvars.state
    fi
fi

# Find the index of the last completed section (looks at compressed velocities)
SECTION=0
for f in out/section*.vel.gz; do
    [ -e "$f" ] || continue
    n=$(basename "$f" | sed -E 's/^section([0-9]+)\.vel(\.gz)?$/\1/')
    if [ "$n" -gt "$SECTION" ]; then
        SECTION=$n
    fi
done

# Refresh out/restart.* from the last completed section if any
if [ "$SECTION" -gt 0 ]; then
    rm -f out/restart.*
    for ext in coor vel xsc xst colvars.state; do
        if [ -f "out/section${SECTION}.${ext}" ]; then
            cp -f "out/section${SECTION}.${ext}" "out/restart.${ext}"
        elif [ -f "out/section${SECTION}.${ext}.gz" ]; then
            gunzip -c "out/section${SECTION}.${ext}.gz" > "out/restart.${ext}"
        fi
    done
fi

# Are we done?
if [ "$SECTION" -ge "$NSECTION" ]; then
    echo "Window $ANALOG/$WIN finished ($SECTION/$NSECTION sections)."
    exit 0
fi

# Run at most SECTIONS_PER_SUBMISSION sections in this single submission
COUNTER=$((SECTION + 1))
LAST_SECTION_THIS_SUBMISSION=$((SECTION + SECTIONS_PER_SUBMISSION))
if [ "$LAST_SECTION_THIS_SUBMISSION" -gt "$NSECTION" ]; then
    LAST_SECTION_THIS_SUBMISSION=$NSECTION
fi

while [ "$COUNTER" -le "$LAST_SECTION_THIS_SUBMISSION" ]; do
    INPUTNAME="out/restart"
    OUTPUTNAME="out/section${COUNTER}"

    # Render the per-section input file from the production template
    sed -e "s|__INPUTNAME__|${INPUTNAME}|g" \
        -e "s|__OUTPUTNAME__|${OUTPUTNAME}|g" \
        production.namd > run_section${COUNTER}.namd

    $NAMD run_section${COUNTER}.namd > out/section${COUNTER}.out

    # Refresh restart files from this section for the next iteration
    rm -f out/restart.*
    for ext in coor vel xsc xst colvars.state; do
        if [ -f "out/section${COUNTER}.${ext}" ]; then
            cp -f "out/section${COUNTER}.${ext}" "out/restart.${ext}"
        fi
    done

    # Compress the bulky binary outputs (mirrors the membrane recipe)
    # Keep .dcd uncompressed for direct VMD/MDAnalysis access during analysis
    # Keep .colvars.traj uncompressed for WHAM postprocessing
    for ext in coor vel xsc xst colvars.state out; do
        if [ -f "out/section${COUNTER}.${ext}" ]; then
            gzip -f "out/section${COUNTER}.${ext}"
        fi
    done

    COUNTER=$((COUNTER + 1))
done

SECTION=$((COUNTER - 1))
if [ "$SECTION" -lt "$NSECTION" ]; then
    echo "Window $ANALOG/$WIN completed $SECTION/$NSECTION sections; submitting next 50-section job."
    sbatch submit.sh
else
    echo "Window $ANALOG/$WIN finished ($SECTION/$NSECTION sections)."
fi
