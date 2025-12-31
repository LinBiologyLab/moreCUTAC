#!/bin/bash
# Copyright (c) 2024 Zheng Wei @ Lin Lab, Fred Hutch Cancer Center
# See LICENSE file in the project root for more information.

if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    if [ -z "$1" ]; then
        echo "SLURM_ARRAY_TASK_ID (job ID) is not set and no parameter provided."
        exit 1
    else
        SLURM_ARRAY_TASK_ID=$1
    fi
fi


workdir=$(awk 'NR==1{print $1}' 00.settings/config.txt)
inputdir=${workdir}/01.SampleList
runID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/runID.txt)
inputdir1=${workdir}/14.genBedGraph
stepdir=${workdir}/15.genBigWig
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.err
if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

module purge
ml Kent_tools/20201201-linux.x86_64

bedGraphToBigWig ${inputdir1}/${runID}.bg  ${refdir}/genome.chrom.sizes  ${stepdir}/${runID}.bw > ${outfile} 2> ${errfile}


bedGraphToBigWig ${inputdir1}/${runID}.norm.bg  ${refdir}/genome.chrom.sizes  ${stepdir}/${runID}.norm.bw > ${outfile} 2> ${errfile}




