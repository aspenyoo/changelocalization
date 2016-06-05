#!/bin/sh
#PBS -o localhost:${PBS_O_WORKDIR}/
#PBS -e localhost:${PBS_O_WORKDIR}/
#PBS -M la67@nyu.edu

NAME=changelocalization

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab/2015a
export MATLABPATH=${MATLABPATH}:/${ROOTPATH}/${NAME}/matlab:${ROOTPATH}/MATLAB
source ${ROOTPATH}/MATLAB/setpath.sh

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
	IID=${PBS_ARRAYID}
fi
#Check if running as an array job
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi

# Run the program

#PARAMS=$(awk "NR==${IID} {print;exit}" ${INPUTFILE})

echo proc ${IID} max run ${RUNMAX} samples [${NSAMPLES},${NSAMPLESFINAL}]

cat<<EOF | matlab -nodisplay -singleCompThread
addpath(genpath('${ROOTPATH}/MATLAB'));
addpath(genpath('${ROOTPATH}/Aspen/${NAME}'));
cd('${WORKDIR}');
runmax=${RUNMAX};
nSamples=[${NSAMPLES},${NSAMPLESFINAL}];
fit_grant_jun2016(${IID},runmax,nSamples);
EOF
