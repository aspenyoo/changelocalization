#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=16:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=4GB
#PBS -m abe
#PBS -N REM

index=${PBS_ARRAYID}
job=${PBS_JOBID}
walltime_lim=${PBS_WALLTIME}
script_name=${PBS_JOBNAME}
module purge
module load matlab

export MATLABPATH=/home/ay963/matlab-scripts
cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/mat-lab-scripts'))
addpath(genpath('/home/ay963/changelocalization’))

subjVec = {‘ALM’,’DR’,’MR’,’EN’};
exptypeVec = {‘Delay’,’Contrast’};

blah = num2str($index);
isubj = str2double(blah(1));
subjid = subjVec{isubj};
model = str2double(blah(2));
exptype = exptypeVec{str2double(blah(3))};


fit_grant_jun2016(subjid,exptype,model)

EOF


