#!/usr/bin/env bash
#
# setup_all.sh
#
# Builds 3 analogs x 38 umbrella windows = 114 simulation directories
# under  umbrella/runs/<analog>/win<NN>/   ,  each containing:
#       system.psf  system.pdb  ref.pdb
#       step6.1_equilibration.namd .. step6.6_equilibration.namd
#       production.namd  colvars.in  submit.sh
#
# Window centers: solute 1 from z = 0 to 37 A, solute 2 always 3.7 nm below.
# (MacCallum 2008, transposed to a single bilayer leaflet.)
#
# Run on the cluster (where VMD is available):
#       bash setup_all.sh
#
# Then launch all jobs with:
#       bash launch_all.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# repository layout (everything below is relative to the project root)
# ---------------------------------------------------------------------------
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

ANALOG_DIR="$HERE/build/analog"
MEMB_DIR="$HERE/build/membrane/charmm-gui/namd"
TPL_DIR="$HERE/templates"
BUILD_TCL="$HERE/build/build_window.tcl"
RUNS_DIR="$HERE/runs"
TOPPAR_DIR="$HERE/toppar"

# ---------------------------------------------------------------------------
# 1. assemble a single toppar/ directory pointed to by every NAMD config
# ---------------------------------------------------------------------------
mkdir -p "$TOPPAR_DIR"
# membrane / water / lipid CHARMM36
#for f in par_all36m_prot.prm par_all36_na.prm par_all36_carb.prm \
#         par_all36_lipid.prm par_all36_cgenff.prm par_interface.prm \
#         toppar_water_ions.str toppar_all36_lipid_cholesterol.str \
#         toppar_all36_lipid_bacterial.str ; do
#    ln -sf "$MEMB_DIR/toppar/$f" "$TOPPAR_DIR/$f"
#done
# analog (side chain) topology + parameters
#ln -sf "$ANALOG_DIR/toppar-namd/top_all36_sidechains.str" \
#       "$TOPPAR_DIR/top_all36_sidechains.str"

# ---------------------------------------------------------------------------
# 2. window definitions
# ---------------------------------------------------------------------------
ANALOGS=("SCL" "SCC" "SCY" "SCS" "SCR" "SCD")

NWIN=38
DZ_PAIR=-37.0         # solute 2 is 3.7 nm BELOW solute 1 along z
Z_MIN=0.0             # first window: solute 1 at z = 0
DZ_STEP=1.0           # 1 A between consecutive windows
                      # window 0  : solute1 z=0,  solute2 z=-37
                      # window 37 : solute1 z=37, solute2 z=0

# ---------------------------------------------------------------------------
# 3. main loop -- build PSF/PDB + render templates for every window
# ---------------------------------------------------------------------------
for ANALOG in "${ANALOGS[@]}"; do
    echo $ANALOG
    al=$(echo "$ANALOG" | tr '[:upper:]' '[:lower:]')
    for ((i=0; i<$NWIN; i++)); do
        WIN=$(printf "%02d" "$i")
        Z1=$(awk "BEGIN{printf \"%.2f\", $Z_MIN + $i * $DZ_STEP}")
        Z2=$(awk "BEGIN{printf \"%.2f\", $Z1 + $DZ_PAIR}")
        DIR="$RUNS_DIR/$al/win$WIN"
        mkdir -p "$DIR"
        echo "==> $al  win=$WIN  z1=$Z1  z2=$Z2"

        # ----- build the system (psfgen + water removal) -------------------
        if [ ! -f "$DIR/system.psf" ]; then
            ( cd "$DIR" && \
              vmd -dispdev text -e "$BUILD_TCL" -args \
                  "$ANALOG" "$Z1" "$DIR" "$MEMB_DIR" "$ANALOG_DIR" \
                  > build.log 2>&1 )
        fi

        # ----- copy / render templates -------------------------------------
        for s in 1 2 3 4 5 6; do
            cp -f "$TPL_DIR/step6.${s}_equilibration.namd" \
                  "$DIR/step6.${s}_equilibration.namd"
        done

        sed -e "s|__WIN__|$WIN|g"   \
            -e "s|__Z1__|$Z1|g"     \
            -e "s|__Z2__|$Z2|g"     \
            "$TPL_DIR/production.namd" > "$DIR/production.namd"

        sed -e "s|__WIN__|$WIN|g"   \
            -e "s|__Z1__|$Z1|g"     \
            -e "s|__Z2__|$Z2|g"     \
            "$TPL_DIR/colvars.in" > "$DIR/colvars.in"

        sed -e "s|__WIN__|$WIN|g"      \
            -e "s|__ANALOG__|$al|g"    \
            "$TPL_DIR/submit.sh" > "$DIR/submit.sh"
        chmod +x "$DIR/submit.sh"
    done
done

echo
echo "Setup complete. ${#ANALOGS[@]} analogs x $NWIN windows = $((${#ANALOGS[@]}*NWIN)) windows."
echo "Launch everything with:   bash $HERE/launch_all.sh"
