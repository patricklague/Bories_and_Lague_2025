#!/usr/bin/env bash

#FAITS: "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: 
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"
#"NONE"

aafile=("SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM")

traj=(1 2 3)


for aa in "${aafile[@]}"
do
  for t in "${traj[@]}"
  do
    cp /media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa/analyses/traj$t/data/thickness.dat ./thickness$t.dat
    if [ $t = 4 ]; then
      mv thickness$t.dat thickness1.dat
    fi
  done
  cat get_thickness.py | sed s=FILENAME=${aa,,}= > get_temp.py
  python get_temp.py
  rm thickness*.dat
  if [ $aa = 'NONE' ]; then
    mv none-thickness.dat popc-thickness.dat
  fi
done

