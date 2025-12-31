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
stepdir=${workdir}/10.callpeakByStatus
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err
if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

module purge
ml MACS2/2.2.9.1-foss-2022b


macs2 callpeak -t ${inputdir1}/${statusID}.bam --format BAM --keep-dup all -q 0.1 --call-summits --name ${statusID} --outdir ${stepdir} -g hs --nomodel --shift -75 --extsize 150  -B --SPMR 2> ${stepdir}/${statusID}_narrow.summary.txt > ${outfile}.narrow

#macs2 callpeak -t ${inputdir1}/${statusID}.bam   --format BAMPE --keep-dup all --call-summits --outdir ${stepdir} -g hs --name ${statusID}  -q 0.1 2> ${stepdir}/${statusID}.summary.txt > ${outfile} 

macs2 callpeak -t  ${inputdir1}/${statusID}.bam -f BAMPE --keep-dup all -n ${statusID} --outdir ${stepdir} -q 0.1 --nolambda --broad  -g hs 2> ${stepdir}/${statusID}.braod.summary.txt > ${outfile}.broad


module purge
ml BEDTools/2.30.0-GCC-11.2.0
bedtools groupby  -g 1,2,3 -c 9 -o max -i ${stepdir}/${statusID}_peaks.narrowPeak |   sort -k4,4nr  |  awk '($1 !~ /_/ && $1 !~ /chrM/) { print }' | bedtools intersect  -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_peaks.narrowPeak.sorted.bed


bedtools groupby  -g 1,2,3 -c 9 -o max -i ${stepdir}/${statusID}_peaks.broadPeak |  sort -k4,4nr   | awk '($1 !~ /_/ && $1 !~ /chrM/) { print }' | bedtools intersect  -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_peaks.broadPeak.sorted.bed
