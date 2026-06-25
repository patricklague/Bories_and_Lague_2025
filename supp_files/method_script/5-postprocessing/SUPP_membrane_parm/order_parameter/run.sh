#!/usr/bin/env bash

aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

chain=(2 3)
traj=(1 2 3)

# Canonical destination for the per-analog merged file
OUTDIR="../../../../../figures/data/SUPP_membrane_parm/order_parameter"
mkdir -p "$OUTDIR"

for i in "${chain[@]}"
do
  for aa in "${aafile[@]}"
  do
    # Original (local) source of the per-trajectory data:
    DIR=../../../POPC-aa/POPC-$aa

    for t in "${traj[@]}"
    do
      cp $DIR/analyses/traj$t/data/orderParameters/orderparameters-chain${i}*.dat .
      for j in 401-600 601-800 801-1000
      do
        mv orderparameters-chain${i}-$j.dat scd${i}-t${t}-$j.dat
      done
    done
    python get_scd.py

    if [ "$aa" = "NONE" ]; then
      out="popc-chain$i.dat"
    else
      out="${aa,,}-chain$i.dat"
    fi
    mv trajectory_scd.dat "$out"
    mv "$out" "$OUTDIR/$out"
    rm scd*-t*.dat
  done
done

