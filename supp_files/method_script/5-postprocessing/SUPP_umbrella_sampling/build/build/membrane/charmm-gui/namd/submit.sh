#! /bin/bash
#SBATCH --account=def-plague
#SBATCH --job-name=BOR
#SBATCH --gpus-per-node=1          # Number of GPU(s) per node
#SBATCH --cpus-per-task=8         # CPU cores/threads
#SBATCH --mem 2048               # memory per node
#SBATCH --time=0-12:00            # time (DD-HH:MM)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# modules may vary according to the stack, check with "module avail"
#module load nixpkgs/16.09  
#module load gcc/5.4.0
#module load cuda/8.0.44
module load namd-multicore/2.14
# Number of DCD sections 
# time of trajectory = $nstep * $nsection * timestep
# nsteps is now set in step5 file
NSECTION=200

# namd command, adjust +p10 to number of available cores as requested to queue system
#NAMD="/home/plague/projects/def-plague/bin/NAMD_Git-2018-06-14_Linux-x86_64-multicore-CUDA/namd2 +p5 +idlepoll"
NAMD="namd2 +p$SLURM_CPUS_PER_TASK +idlepoll"

#
# Do not edit below...
#

# need to start the trajectory?
if [ ! -f out/section1.dcd ]; then

    # start trajectory
    $NAMD step6.1_equilibration.inp >& step6.1_equilibration.out
    $NAMD step6.2_equilibration.inp >& step6.2_equilibration.out
    $NAMD step6.3_equilibration.inp >& step6.3_equilibration.out
    $NAMD step6.4_equilibration.inp >& step6.4_equilibration.out
    $NAMD step6.5_equilibration.inp >& step6.5_equilibration.out
    $NAMD step6.6_equilibration.inp >& step6.6_equilibration.out

    mkdir -p out
    cp step6.6_equilibration.xst out/restart.xst
    cp step6.6_equilibration.xsc out/restart.xsc
    cp step6.6_equilibration.coor out/restart.coor
    cp step6.6_equilibration.vel out/restart.vel

fi

#
# restart trajectory
#

# find the number of the last section
SECTION=0
if [ -f out/section1.dcd ]; then
    for file in `ls out/section*out.gz`;
        do
        temp=`echo $file | awk -F . '{ print $1 }' | awk -F on '{ print $2 }'`
        if  [ "$temp" -gt "$SECTION" ]; then
            SECTION=$temp
        fi
    done
    echo Section $SECTION will be used to restart the trajectory
fi

# make sure we have the right files for restart
if [ $SECTION -gt 0 ]; then

    rm -f out/restart*
    cp -f out/section$SECTION.xst.gz out/restart.xst.gz
    cp -f out/section$SECTION.coor.gz out/restart.coor.gz
    cp -f out/section$SECTION.vel.gz out/restart.vel.gz
    cp -f out/section$SECTION.xsc.gz out/restart.xsc.gz
    gunzip out/restart*
fi
let COUNTER=SECTION+1

# create appropriate input file using step5_production.inp as template
inputname="out/restart"
outputname="out/section$COUNTER"
sed "s/step6.6_equilibration/out\/restart/" step7_production.inp | \
  sed "s/step7_production/out\/section$COUNTER/" > step7_run.inp

# run the simulation for 1 nanosecond
$NAMD step7_run.inp > out/section$COUNTER.out
gzip out/section$COUNTER.xst
gzip out/section$COUNTER.coor
#gzip out/section$COUNTER.dcd
gzip out/section$COUNTER.vel
gzip out/section$COUNTER.xsc
gzip out/section$COUNTER.out

if [ $COUNTER -le $NSECTION ]; then
  /opt/software/slurm/bin/sbatch submit.sh > jobid
fi

exit

