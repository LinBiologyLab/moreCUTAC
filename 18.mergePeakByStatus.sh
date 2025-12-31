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
inputdir1=${workdir}/07.callpeak
inputdir2=${workdir}/10.callpeakByStatus
stepdir=${workdir}/18.mergePeakByStatus
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
ml BEDTools/2.30.0-GCC-11.2.0

cat ${inputdir}/runID.txt | grep ${statusID} | awk -v inputdir1=$inputdir1  '{print inputdir1"/"$1"_pe_peaks.narrowPeak.sorted"}'  | xargs -n 1  cat  | sort -k1,1 -k2,3n |  awk '($1 !~ /_/ && $1 !~ /chrM/) { print }'  | bedtools merge -i stdin | bedtools intersect -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_pe_peaks.narrowPeak.merged.bed

cat ${inputdir}/runID.txt | grep ${statusID} | awk -v inputdir1=$inputdir1  '{print inputdir1"/"$1"_pe_peaks.narrowPeak.sorted"}'  | xargs -n 1   head -n 50000 | sort -k1,1 -k2,3n | awk '($1 !~ /_/ && $1 !~ /chrM/) { print }'  | bedtools merge -i stdin | bedtools intersect -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_pe_top50kpeaks.narrowPeak.merged.bed

cat ${inputdir}/runID.txt | grep ${statusID} | awk -v inputdir1=$inputdir1  '{print inputdir1"/"$1"_se_peaks.narrowPeak.sorted"}'  | xargs -n 1  cat  | sort -k1,1 -k2,3n |  awk '($1 !~ /_/ && $1 !~ /chrM/) { print }'  | bedtools merge -i stdin | bedtools intersect -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_se_peaks.narrowPeak.merged.bed

cat ${inputdir}/runID.txt | grep ${statusID} | awk -v inputdir1=$inputdir1  '{print inputdir1"/"$1"_se_peaks.narrowPeak.sorted"}'  | xargs -n 1   head -n 50000 | sort -k1,1 -k2,3n | awk '($1 !~ /_/ && $1 !~ /chrM/) { print }'  | bedtools merge -i stdin | bedtools intersect -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_se_top50kpeaks.narrowPeak.merged.bed


cat ${inputdir}/runID.txt | grep ${statusID} | awk -v inputdir1=$inputdir1  '{print inputdir1"/"$1"_peaks.broadPeak.sorted"}'  | xargs -n 1 cat | sort -k1,1 -k2,3n | awk '($1 !~ /_/ && $1 !~ /chrM/) { print }' | bedtools merge -i stdin | bedtools intersect -a stdin -b ${refdir}/blacklist.rev.bed > ${stepdir}/${statusID}_peaks.broadPeak.merged.bed

#cat ${inputdir2}/${statusID}_peaks.narrowPeak.sorted.bed  ${stepdir}/${statusID}_pe_peaks.narrowPeak.merged.bed ${stepdir}/${statusID}_se_peaks.narrowPeak.merged.bed  ${stepdir}/${statusID}_peaks.broadPeak.merged.bed | awk 'OFS="\t"{print $1,$2,$3}' | sort -k1,1 -k2,3n | bedtools merge -i stdin > ${stepdir}/${statusID}_peaks.merged.bed

cat ${inputdir2}/${statusID}_peaks.narrowPeak.sorted.bed  ${inputdir2}/${statusID}_peaks.broadPeak.sorted.bed  | awk 'OFS="\t"{print $1,$2,$3}' | sort -k1,1 -k2,3n | bedtools merge -i stdin > ${stepdir}/${statusID}_peaks.merged.bed
