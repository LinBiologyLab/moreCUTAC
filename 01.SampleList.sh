#!/bin/bash
# Copyright (c) 2024 Zheng Wei @ Lin Lab, Fred Hutch Cancer Center
# See LICENSE file in the project root for more information.
workdir=$(awk 'NR==1{print $1}' 00.settings/config.txt)
stepdir=${workdir}/01.SampleList
mkdir -p ${stepdir}
datadir=$(awk 'NR==2{print $1}' 00.settings/config.txt)

cp -r 00.settings ${workdir}
echo "Do not make changes to files In this folder. If you want to rerun the script with updated settings and input datalist, please go the 00.settings in correspoinding script folder" > ${workdir}/00.settings/DoNotMakeChangesToFilesInThisFolder.txt

metadatafile=00.settings/metadata.txt

cat $metadatafile | awk 'NR > 1'  > ${stepdir}/datalist.txt
awk -F'\t' '{print $4"_"$5"_"$6"_"$7}' ${stepdir}/datalist.txt | sort -k1,1 | uniq > ${stepdir}/repID.txt
awk -F'\t' '{print $4"_"$5"_"$6}' ${stepdir}/datalist.txt | sort -k1,1 | uniq > ${stepdir}/statusID.txt
awk -F'\t' '{print $4"_"$5}' ${stepdir}/datalist.txt | sort -k1,1 | uniq > ${stepdir}/sampleID.txt
awk -F'\t' '{print $4"_"$5"_"$6"_"$7"_"$8"_"$9}' ${stepdir}/datalist.txt  | sort -k1,1  > ${stepdir}/runID.txt
awk -F'\t' 'OFS="\t"{print $4"_"$5"_"$6"_"$7"_"$8"_"$9,$1,$2,$11}' ${stepdir}/datalist.txt | sort -k1,1 |  awk -F'+' 'OFS="\t"{print $1,$2}' > ${stepdir}/fastqrun.txt
