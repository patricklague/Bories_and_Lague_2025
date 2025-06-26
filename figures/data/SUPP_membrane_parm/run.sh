#!/usr/bin/env bash

#FAITS: "SCV" "GLYD" "SCA" "SCP" "SCW" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC"
#A FAIRE: "SCM"
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"

aafile=( "SCM")

traj=(1 2 3)


for aa in "${aafile[@]}"
do
  mkdir popc/distributions/${aa,,}
  for t in "${traj[@]}"
  do
    cp ../results/homoPOPC-aa/homoPOPC-$aa/analyses/traj$t/data/densityProfiles/profile-aa-600.dat ./trajectory$t.dat
    if [ $t = 4 ]; then
      mv trajectory$t.dat trajectory1.dat
    fi
  done
  python pmf-from-distribution5.py
  mv trajectory*dat popc/distributions/${aa,,}/
  mv pmf_moyen.dat popc/pmfs/${aa,,}.dat
  mv pmf.dat pmf-data/pmf_${aa,,}.dat
done

