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
inputdir1=${workdir}/03.mapping
stepdir=${workdir}/04.addRG
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${runID}.err

if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi


tecrepID=$(echo ${runID} | awk -F'_' '{print $1"_"$2"_"$3"_"$4"_"$5}' )
statusID=$(echo ${runID} | awk -F'_' '{print $1"_"$2"_"$3}' )
bindir=$(awk 'NR==4{print $1}' 00.settings/config.txt)

module purge
#ml Java/17.0.6

#java -jar ${bindir}/picard/build/libs/picard.jar AddOrReplaceReadGroups -I ${inputdir1}/${runID}.bam    -O ${stepdir}/${runID}.bam     -SORT_ORDER coordinate    -RGID ${runID}     -RGLB ${tecrepID}     -RGPL ILLUMINA    -RGSM ${statusID}   -RGPU ${runID}     --CREATE_INDEX true > ${outfile} 2> ${errfile}

ml picard/2.25.1-Java-11
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I ${inputdir1}/${runID}.bam    -O ${stepdir}/${runID}.bam     -SORT_ORDER coordinate    -RGID ${runID}     -RGLB ${tecrepID}     -RGPL ILLUMINA    -RGSM ${statusID}   -RGPU ${runID}     --CREATE_INDEX true > ${outfile} 2> ${errfile}
