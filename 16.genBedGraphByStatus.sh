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
statusID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/statusID.txt)
inputdir1=${workdir}/09.removealldupByStatus
stepdir=${workdir}/16.genBedGraphByStatus
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err
if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

statusID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/statusID.txt)

module purge
ml BEDTools/2.30.0-GCC-11.2.0

bedtools genomecov -bg  -ibam ${inputdir1}/${statusID}.bam | sort -k1,1 -k2,3n  > ${stepdir}/${statusID}.bg  2> ${errfile}

allcounts=$(cat ${inputdir1}/${statusID}.count.txt | awk '{print $1/1000000}')

bedtools genomecov -bg  -ibam ${inputdir1}/${statusID}.bam | sort -k1,1 -k2,3n | awk -v allcounts=${allcounts} 'OFS="\t"{print $1, $2, $3, $4/allcounts}' > ${stepdir}/${statusID}.norm.bg  2>> ${errfile}


bedtools merge -i ${stepdir}/${statusID}.bg > ${stepdir}/${statusID}.cov_reads.bed



