#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=16:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=4GB
#PBS -m abe
#PBS -N VSTM_CL

index=${PBS_ARRAYID}
job=${PBS_JOBID}
walltime_lim=${PBS_WALLTIME}s
script_name=${PBS_JOBNAME}
module purge
module load matlab

export MATLABPATH=/home/ay963/matlab-scripts
cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/mat-lab-scripts'))
addpath(genpath('/home/ay963/changelocalization’))


blah = num2str($index);
model = str2double(blah(1));
isubj = str2double(blah(2));
runlist = str2double(blah(3:end));

fit_grant_jun2016(isubj,model,runlist)

EOF


