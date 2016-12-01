#!/bin/bash
PROJECT=changelocalization
SHORTNAME=BCL
ROOTPATH="/home/la67"
BASEDIR="${ROOTPATH}/Aspen/${PROJECT}"
SOURCEDIR="${BASEDIR}/matlab"
JOBSCRIPT="${BASEDIR}/scripts/fitgrant.sh"

#Job parameters
RUN=${1}
WORKDIR="/scratch/la67/${PROJECT}/run${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}
MAXID=400
RUNTIME=12:00:00
MAXRT=NaN
RUNMAX=50
NSAMPLES=20000
NSAMPLESFINAL=500000

#Job list is second argument
if [[ ! -z "$2" ]]; then
        JOBLIST=$2
else
	JOBLIST="1-${MAXID}"
fi

#RESOURCES="nodes=1:ppn=1,mem=4GB,walltime=${RUNTIME},feature=ivybridge_20p_64GB_3000MHz"
RESOURCES="nodes=1:ppn=1,mem=8GB,walltime=${RUNTIME}"

#Convert from spaces to commas
JOBLIST=${JOBLIST// /,}
echo JOBS $JOBLIST

JOBNAME=${SHORTNAME}${RUN}
qsub -t ${JOBLIST} -v PROJECT=${PROJECT},ROOTPATH=${ROOTPATH},MAXID=$MAXID,MAXRT=$MAXRT,WORKDIR=$WORKDIR,RUNMAX=${RUNMAX},NSAMPLES=${NSAMPLES},NSAMPLESFINAL=${NSAMPLESFINAL} -l ${RESOURCES} -N ${JOBNAME} ${JOBSCRIPT}
