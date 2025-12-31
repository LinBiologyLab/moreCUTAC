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
statusID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/statusID.txt)
inputdir1=${workdir}/19.calibration
inputdir2=${workdir}/18.mergePeakByStatus
stepdir=${workdir}/20.mutect2TumorOnly
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err

if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi



runbam=$(grep ${statusID} ${inputdir}/runID.txt  | awk -v inputdir1=${inputdir1} '{print "-I "inputdir1"/"$1".bam"}')


module purge
ml GATK/4.4.0.0-GCCcore-12.2.0-Java-17


refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

gatk Mutect2  -R ${refdir}/genome.fa ${runbam}  -tumor   ${statusID}  -pon ${refdir}/gatk/pon.vcf.gz   --germline-resource ${refdir}/gatk/af-only-gnomad.vcf.gz   -L ${inputdir2}/${statusID}_peaks.merged.bed --f1r2-tar-gz ${stepdir}/${statusID}.f1r2.tar.gz    -O ${stepdir}/${statusID}.mutect2.vcf    -bamout ${stepdir}/${statusID}.mutect2.bam --native-pair-hmm-threads ${cpucores} --tmp-dir ${stepdir}   > ${outfile} 2> ${errfile}

module purge
ml SAMtools/1.14-GCC-11.2.0

samtools index ${stepdir}/${statusID}.mutect2.bam

