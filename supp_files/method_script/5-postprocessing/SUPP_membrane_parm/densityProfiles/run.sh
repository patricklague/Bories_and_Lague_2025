#!/usr/bin/env bash

aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

profiles=("chains" "choline" "phosphate" "total" "water")
traj=(1 2 3)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../../figures/data/SUPP_membrane_parm/densityProfiles"
mkdir -p "$OUTDIR"

for i in "${profiles[@]}"
do
  for aa in "${aafile[@]}"
  do
    # Original (local) source of the per-trajectory data:
    DIR=../../../POPC-aa/POPC-$aa

    for t in "${traj[@]}"
    do
      cp $DIR/analyses/traj$t/data/densityProfiles/profile-${i}*.dat .
      for j in 401-600 601-800 801-1000
      do
        mv profile-${i}-$j.dat trajectory$t-$j.dat
      done
    done
    python get_density.py

    out="${aa,,}-$i.dat"
    if [ "$aa" = "NONE" ]; then
      out="popc-$i.dat"
    fi
    mv trajectory-dens.dat "$out"
    mv "$out" "$OUTDIR/$out"
    rm trajectory*.dat
  done
done

# When all densities have been calculated
#for i in "${profiles[@]}"
#do
#  cat calculate_mean_density.py | sed s=PROFILE=$i= > temp.py
#  python temp.py
#  mv all_dens.dat density-$i.dat
#done
