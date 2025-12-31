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
inputdir1=${workdir}/06.removealldup
stepdir=${workdir}/07.callpeak
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.err
if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

module purge
ml MACS2/2.2.9.1-foss-2022b

macs2 callpeak -t ${inputdir1}/${runID}.bam   --format BAMPE --keep-dup all --call-summits --outdir ${stepdir} -g hs --name ${runID}_pe  -q 0.1 2> ${stepdir}/${runID}.pe_narrow.summary.txt > ${outfile}.pe_narrow

sort -k9,9nr ${stepdir}/${runID}_pe_peaks.narrowPeak | awk '($1 !~ /_/ && $1 !~ /chrM/) { print }'  > ${stepdir}/${runID}_pe_peaks.narrowPeak.sorted


macs2 callpeak -t ${inputdir1}/${runID}.bam --format BAM --keep-dup all -q 0.1 --call-summits --name ${runID}_se --outdir ${stepdir} -g hs --nomodel --shift -75 --extsize 150  -B --SPMR 2> ${stepdir}/${runID}_se_narrow.summary.txt > ${outfile}.se_narrow

sort -k9,9nr ${stepdir}/${runID}_se_peaks.narrowPeak |  awk '($1 !~ /_/ && $1 !~ /chrM/) { print }' > ${stepdir}/${runID}_se_peaks.narrowPeak.sorted


macs2 callpeak -t  ${inputdir1}/${runID}.bam -f BAMPE --keep-dup all -n ${runID} --outdir ${stepdir} -q 0.1 --nolambda --broad  -g hs 2> ${stepdir}/${runID}.braod.summary.txt > ${outfile}.broad 

sort -k9,9nr ${stepdir}/${runID}_peaks.broadPeak | awk '($1 !~ /_/ && $1 !~ /chrM/) { print }' > ${stepdir}/${runID}_peaks.broadPeak.sorted






