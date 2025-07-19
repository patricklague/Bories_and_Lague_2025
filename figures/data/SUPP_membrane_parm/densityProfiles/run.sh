#!/usr/bin/env bash

#FAITS: "SCYM" "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM"
#A FAIRE: 
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE"
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW"

aafile=( "SCM" )
profiles=("chains" "choline" "phosphate" "total" "water")
traj=(1 2 3)

for i in "${profiles[@]}"
do
  for aa in "${aafile[@]}"
  do
    for t in "${traj[@]}"
    do
      cp /media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa/analyses/traj$t/data/densityProfiles/profile-${i}600.dat ./trajectory$t.dat
      if [ $t = 4 ]; then
        mv trajectory$t.dat trajectory1.dat
      fi
    done
    python get_density.py
    mv trajectory-dens.dat ${aa,,}-$i.dat
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
