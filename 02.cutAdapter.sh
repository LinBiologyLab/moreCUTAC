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
stepdir=${workdir}/02.cutAdapter
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.err

if [[ ! -d "${stepdir}" ]]; then
    mkdir ${stepdir}
fi

if [[ ! -d "${stepdir}/out" ]]; then
    mkdir ${stepdir}/out
fi

runID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/runID.txt)
adapter1=$(awk -v  runID=${runID}  '$1 == runID {print $4}' ${inputdir}/fastqrun.txt | grep '.' | awk '{print "-a "$1}')
adapter2=$(awk -v  runID=${runID}  '$1 == runID {print $5}' ${inputdir}/fastqrun.txt | grep '.' | awk '{print "-A "$1}')

datadir=$(awk 'NR==2{print $1}' 00.settings/config.txt)

fastqfile1=$(awk -v  runID=${runID} -v datadir=${datadir}  '$1 == runID {print datadir"/"$2" "}' ${inputdir}/fastqrun.txt)
fastqfile2=$(awk -v  runID=${runID} -v datadir=${datadir}  '$1 == runID {print datadir"/"$3" "}' ${inputdir}/fastqrun.txt)

module purge
ml FastQC/0.11.8-Java-1.8

if [[ ! -d "${stepdir}/fastqc" ]]; then
    mkdir ${stepdir}/fastqc
fi

fastqc -t ${cpucores} -o ${stepdir}/fastqc $fastqfile1 $fastqfile2 > ${outfile}.fastqc  2> ${errfile}.fastqc
fq1=$(basename $fastqfile1 .fastq.gz)
fq2=$(basename $fastqfile2 .fastq.gz)
#cp ${stepdir}/fastqc/${fq1}_fastqc.html  ${stepdir}/fastqc/${runID}_1_fastqc.html   
#cp ${stepdir}/fastqc/${fq2}_fastqc.html  ${stepdir}/fastqc/${runID}_2_fastqc.html 
#cp ${stepdir}/fastqc/${fq1}_fastqc.zip  ${stepdir}/fastqc/${runID}_1_fastqc.zip
#cp ${stepdir}/fastqc/${fq2}_fastqc.zip  ${stepdir}/fastqc/${runID}_2_fastqc.zip
#ln -s ${stepdir}/fastqc/${fq1}_fastqc.html ${runID}_1_fastqc.html 
#ln -s ${stepdir}/fastqc/${fq2}_fastqc.html ${runID}_2_fastqc.html
#ln -s ${stepdir}/fastqc/${fq1}_fastqc.zip ${runID}_1_fastqc.zip 
#ln -s ${stepdir}/fastqc/${fq2}_fastqc.zip ${runID}_2_fastqc.zip 

module purge
ml cutadapt/4.1-GCCcore-11.2.0


echo $fastqfile1
echo ${stepdir}/${runID}_1.fastq.gz

cutadapt -j ${cpucores} --nextseq-trim 20 -m 20 $adapter1 $adapter2 -Z -o ${stepdir}/${runID}_1.fastq.gz -p ${stepdir}/${runID}_2.fastq.gz $fastqfile1 $fastqfile2 > ${stepdir}/${runID}.txt 2> ${errfile}

olddir=$(pwd)
cd ${stepdir}/fastqc

ln -s ${fq1}_fastqc.html ${runID}_1_fastqc.html
ln -s ${fq2}_fastqc.html ${runID}_2_fastqc.html
ln -s ${fq1}_fastqc.zip ${runID}_1_fastqc.zip
ln -s ${fq2}_fastqc.zip ${runID}_2_fastqc.zip

cd $olddir
