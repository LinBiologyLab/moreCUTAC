bash 01.SampleList.sh

workdir=$(awk 'NR==1{print $1}' 00.settings/config.txt)
cpucores=$(awk 'NR==5{print $1}' 00.settings/config.txt)
memory=$(awk 'NR==6{print $1}' 00.settings/config.txt)
inputdir=${workdir}/01.SampleList
sampleCount=$(wc -l ${inputdir}/sampleID.txt| awk '{print $1}')

#echo "[1-${sampleCount}]"


seq 1 ${sampleCount} | xargs -n 1 bash single_sample.sh

