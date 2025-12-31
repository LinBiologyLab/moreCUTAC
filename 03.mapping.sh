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
cpucores=$(awk 'NR==5{print $1}' 00.settings/config.txt)
inputdir=${workdir}/01.SampleList
runID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/runID.txt)
inputdir1=${workdir}/02.cutAdapter
stepdir=${workdir}/03.mapping
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.err

if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

runID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/runID.txt)

refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

module purge
ml BWA/0.7.17-GCCcore-11.2.0
ml SAMtools/1.16.1-GCC-11.2.0

bwa mem -t ${cpucores} ${refdir}/bwa/genome.fa ${inputdir1}/${runID}_1.fastq.gz ${inputdir1}/${runID}_2.fastq.gz 2> ${errfile}  | samtools view -Sb > ${stepdir}/${runID}.bam  2> ${errfile}.samtools

