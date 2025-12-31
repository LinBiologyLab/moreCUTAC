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
statusID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/statusID.txt )
inputdir1=${workdir}/20.mutect2TumorOnly
stepdir=${workdir}/23.filterMutectCalls
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err

if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi



runID=$(grep $statusID ${inputdir}/runID.txt  | awk -v inputdir2=${inputdir2} '{print "--contamination-table  "inputdir2"/"$1".table"}' )
    
module purge
ml GATK/4.4.0.0-GCCcore-12.2.0-Java-17



refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

isFresh=$(echo $runID | grep 'Fresh' | wc -l | awk '$1>0{print "Fresh"}')

if [ "$isFresh" = "Fresh" ]; then
  gatk FilterMutectCalls    -V ${inputdir1}/${statusID}.mutect2.vcf -R ${refdir}/genome.fa  -O ${stepdir}/${statusID}.vcf > ${outfile}.FilterMutectCalls 2> ${errfile}.FilterMutectCalls
else
  gatk LearnReadOrientationModel \
	  -I ${inputdir1}/${statusID}.f1r2.tar.gz \
	    -O ${stepdir}/${statusID}.read-orientation-model.tar.gz > ${outfile}.LearnReadOrientationModel 2> ${errfile}.LearnReadOrientationModel

  gatk FilterMutectCalls    -V ${inputdir1}/${statusID}.mutect2.vcf -R ${refdir}/genome.fa -ob-priors ${stepdir}/${statusID}.read-orientation-model.tar.gz -O ${stepdir}/${statusID}.vcf > ${outfile}.FilterMutectCalls 2> ${errfile}.FilterMutectCalls

fi



