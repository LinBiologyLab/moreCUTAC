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
inputdir1=${workdir}/08.mergeBamByStatus
stepdir=${workdir}/09.removealldupByStatus
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err

if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi


module purge
ml GATK/4.4.0.0-GCCcore-12.2.0-Java-17


gatk  --java-options "-Xmx500G" MarkDuplicatesSpark  -I ${inputdir1}/${statusID}.bam -O ${stepdir}/${statusID}.bam --remove-all-duplicates --create-output-bam-index true  --tmp-dir ${stepdir}  > ${outfile} 2> ${errfile}

module purge
ml SAMtools/1.14-GCC-11.2.0

samtools view -c -f 3  ${stepdir}/${statusID}.bam  > ${stepdir}/${statusID}.count.txt 2>> ${errfile}
