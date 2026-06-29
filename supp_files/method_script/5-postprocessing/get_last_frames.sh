#!/usr/bin/env bash

aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" "SCF-1" "SCY-1" "SCL-1" "SCV-1" "SCI-1" "SCM-1" \
  "NONE" \
)

traj=(1 2 3)


for aa in "${aafile[@]}"
do
  mkdir ../../output_systems/${aa,,}
  for t in "${traj[@]}"
  do
    cp ../../results/POPC-aa/POPC-$aa/charmm-gui/namd/out$t/section1000.dcd ./trajectory$t.dcd
    cp trajectory$t.dcd traj.dcd
    cp ../../results/POPC-aa/POPC-$aa/charmm-gui/namd/step5_input.psf ./input.psf
    vmd -dispdev text -e get_frame.tcl
    mv last_frame.pdb ../../output_systems/${aa,,}/end_traj$t.pdb
  done
  rm trajectory* input* traj.dcd
done

