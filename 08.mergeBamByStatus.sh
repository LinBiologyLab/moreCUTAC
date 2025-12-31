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
inputdir1=${workdir}/05.removedup
stepdir=${workdir}/08.mergeBamByStatus
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err
mkdir ${stepdir}
mkdir ${stepdir}/out



bindir=$(awk 'NR==4{print $1}' 00.settings/config.txt)

module purge
#ml Java/17.0.6

#cat  ${inputdir}/runID.txt | grep ${statusID} | awk -v datadir=${inputdir1} '{print "-I "datadir"/"$1".bam"}' | xargs -n 100 java -jar ${bindir}/picard/build/libs/picard.jar MergeSamFiles -O ${stepdir}/${statusID}.bam  > ${outfile} 2> ${errfile}

ml picard/2.25.1-Java-11
cat  ${inputdir}/runID.txt | grep ${statusID} | awk -v datadir=${inputdir1} '{print "-I "datadir"/"$1".bam"}' | xargs -n 100 java -jar $EBROOTPICARD/picard.jar MergeSamFiles -O ${stepdir}/${statusID}.bam  > ${outfile} 2> ${errfile}


module purge 
ml SAMtools/1.14-GCC-11.2.0
samtools index  ${stepdir}/${statusID}.bam
samtools view -f 3 -c   ${stepdir}/${statusID}.bam  > ${stepdir}/${statusID}.count.txt 2>> ${errfile}

#ml BEDTools/2.30.0-GCC-11.2.0
#samtools view -f 3 -h  ${stepdir}/${statusID}.bam  | samtools view -h  -SB - | samtools sort -n -T ${stepdir}  | bedtools bamtobed -bedpe -mate1 -i stdin | awk 'OFS="\t"{if($9=="+") print $1,$2+4,$6-5,$7,"0","+"; else print $1,$5+4,$3-5,$7,"0","-"}' | grep -v 'chrM' | grep -v '_'  | sort -k1,1 -k2,3n |  awk '$3>$2&&OFS="\t"{print}' > ${stepdir}/${statusID}.bed



