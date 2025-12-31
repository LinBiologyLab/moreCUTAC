step=$(basename $1 .sh)     # the step bash: xxx.sh 
workdir=$(awk 'NR==1{print $1}' 00.settings/config.txt)
jobcount=$(grep $1  steplist.txt | awk -v workdir=$workdir '{print workdir"/01.SampleList/"$2"ID.txt"}' | xargs -n 1 wc -l | awk '{print $1}')
#jobcount=$2                 # the number of RUN/STATUS  01
cpu=$(awk 'NR==5{print $1}' 00.settings/config.txt)                     # CPU job
mem=$(awk 'NR==6{print $1}' 00.settings/config.txt)                      # Memory per job

sbatch --array=[1-${jobcount}] --nodes=1  --output=out/${step}/%a-%A-%J.out --error=out/${step}/%a-%A-%J.err --mem=${mem}000 --cpus-per-task=${cpu} ${1}
