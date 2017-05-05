#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=16:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=4GB
#PBS -m abe
#PBS -N vstm_cl

index=${PBS_ARRAYID}
job=${PBS_JOBID}
walltime_lim=${PBS_WALLTIME}
script_name=${PBS_JOBNAME}
module purge
module load matlab

export MATLABPATH=/home/ay963/matlab-scripts
cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/changelocalization'))

blah = num2str($index);
model = str2double(blah(1));
isubj = blah(2);
runlist = str2double(blah(3:end));

fit_grant_jun2016(isubj,model,runlist)

EOF


