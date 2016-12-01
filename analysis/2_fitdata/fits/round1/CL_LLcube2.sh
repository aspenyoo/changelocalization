#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=8GB
#PBS -m abe
#PBS -N sep_CL_LLmat2

index=${PBS_ARRAYID}
job=${PBS_JOBID}
walltime_lim=${PBS_WALLTIME}
script_name=${PBS_JOBNAME}
module purge
module load matlab

export MATLABPATH=/home/ay963/matlab-scripts
cat<<EOF | matlab -nodisplay
addpath('/home/ay963/job-scripts')
addpath(genpath('/home/ay963/matlab-scripts'))

subjids = {'PM','SM','DG','SC','AD'}
nSamples = 1e5;

% get model and subject idx
blah = num2str($index);
model = str2num(blah(1));
subjidx = str2num(blah(2));
subjid = subjids{subjidx};

% min 0.5 max 30
gridMat = [-0.69 3.5 30; -0.69 3.5 30; -0.69 3.5 30];
parts = [50 str2num(blah(3:end))];

LLcube(subjid, model, gridMat, nSamples, parts)

EOF


