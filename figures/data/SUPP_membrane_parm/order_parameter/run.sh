#!/usr/bin/env bash

#FAITS: "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: "SCM"
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"

aafile=("NONE" "SCM")
chain=(2 3)
traj=(1 2 3)

for i in "${chain[@]}"
do
  for aa in "${aafile[@]}"
  do
    for t in "${traj[@]}"
    do
      cp /media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa/analyses/traj$t/data/orderParameters/orderparameters-chain${i}600.dat ./trajectory$t.dat
      if [ $t = 4 ]; then
        mv trajectory$t.dat trajectory1.dat
      fi
    done
    python get_scd.py
    if [ $aa = "NONE" ]; then
      mv trajectory_scd.dat popc-chain$i.dat
    else
      mv trajectory_scd.dat ${aa,,}-chain$i.dat
    fi
    rm trajectory*.dat
  done
done

