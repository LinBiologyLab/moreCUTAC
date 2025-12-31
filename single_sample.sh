#!/bin/sh

workdir=$(awk 'NR==1{print $1}' 00.settings/config.txt)
cpucores=$(awk 'NR==5{print $1}' 00.settings/config.txt)
memory=$(awk 'NR==6{print $1}' 00.settings/config.txt)
inputdir=${workdir}/01.SampleList
sampleID=$(awk -v  bbb=${SLURM_ARRAY_TASK_ID} 'FNR == bbb {print}' ${inputdir}/sampleID.txt)

grep -nr ${sampleID} ${inputdir}/runID.txt  | awk -F':' '{print $1}' | xargs -n 1 -P 5 bash script_by_run.sh
grep -nr ${sampleID} ${inputdir}/statusID.txt | awk -F':' '{print $1}' | xargs -n 1 -P 5 bash script_by_status.sh

