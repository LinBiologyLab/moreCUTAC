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
inputdir1=${workdir}/23.filterMutectCalls
stepdir=${workdir}/24.selectVariants
mkdir ${stepdir}
mkdir ${stepdir}/out
outfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.out
errfile=${stepdir}/out/${SLURM_ARRAY_TASK_ID}.${statusID}.err


sampleID=$(echo $statusID | awk -F'_' '{print $1"_"$2}')
    
module purge
ml GATK/4.4.0.0-GCCcore-12.2.0-Java-17



refdir=$(awk 'NR==3{print $1}' 00.settings/config.txt)

isFresh=$(echo $statusID | grep 'Fresh' | wc -l | awk '$1>0{print "Fresh"}')

if [ "$isFresh" = "Fresh" ]; then
gatk SelectVariants \
	 -R ${refdir}/genome.fa \
	  -V ${inputdir1}/${statusID}.vcf \
	   --exclude-filtered true \
	    -O ${stepdir}/${statusID}.mutect2.vcf > ${outfile} 2> ${errfile}
else
  

gatk SelectVariants \
	 -R ${refdir}/genome.fa \
	  -V ${inputdir1}/${statusID}.vcf \
	   --exclude-filtered true \
	    -O ${stepdir}/${statusID}.vcf > ${outfile} 2> ${errfile}

grep '#' ${stepdir}/${statusID}.vcf > ${stepdir}/${statusID}.mutect2.vcf

grep -v '#' ${stepdir}/${statusID}.vcf > ${stepdir}/${statusID}.mutect2.tmp
awk '{print $10}'  ${stepdir}/${statusID}.mutect2.tmp  | awk -F':' '{print $5","$6}'  | awk -F',' '{print $2+$4}'  |  paste  ${stepdir}/${statusID}.mutect2.tmp - | awk -F'\t' '!((($4~"^C$"&&$5~"^A|T$")||($4~"^G$"&&$5~"^A|T$"))&&$11~"^1$")&&OFS="\t"{print}' >> ${stepdir}/${statusID}.mutect2.vcf


fi

module purge
ml BEDTools/2.30.0-GCC-11.2.0

if [[ "$statusID" == *"_D"* ]]; then
    cat ${stepdir}/${statusID}.mutect2.vcf | grep -v 'chrM' |  grep -v 'chrY' | grep -v '_alt' | grep -v '_random' | grep -v '_fix' | grep -v 'Un_'  | bedtools intersect -v -a stdin -b ${refdir}/blacklist.bed  | wc -l | awk -v sampleID=${sampleID} 'OFS="\t"{print $1,sampleID,"Tonly"}' >  ${stepdir}/${statusID}.mutect2.txt
fi


if [[ "$statusID" == *"_N"* ]]; then
    cat ${stepdir}/${statusID}.mutect2.vcf | grep -v 'chrM' |  grep -v 'chrY' | grep -v '_alt' | grep -v '_random' | grep -v '_fix' | grep -v 'Un_'  | bedtools intersect -v -a stdin -b ${refdir}/blacklist.bed  | wc -l | awk -v sampleID=${sampleID} 'OFS="\t"{print $1,sampleID,"Nonly"}' >  ${stepdir}/${statusID}.mutect2.txt
fi

