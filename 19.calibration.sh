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
inputdir1=${workdir}/05.removedup
stepdir=${workdir}/19.calibration
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.err
if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

module purge
ml GATK/4.4.0.0-GCCcore-12.2.0-Java-17

refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

gatk --java-options "-Xmx100G"  BaseRecalibrator  -I ${inputdir1}/${runID}.bam -R ${refdir}/genome.fa  -O ${stepdir}/${runID}.table --known-sites ${refdir}/gatk/af-only-gnomad.vcf.gz  > ${outfile} 2> ${errfile}
 
gatk --java-options "-Xmx100G" ApplyBQSR -R  ${refdir}/genome.fa  -I ${inputdir1}/${runID}.bam --bqsr-recal-file  ${stepdir}/${runID}.table -O ${stepdir}/${runID}.bam  > ${outfile} 2> ${errfile}

