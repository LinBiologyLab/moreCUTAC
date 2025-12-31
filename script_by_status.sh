statusIDNo=$1
steplist=$(grep -v '#' steplist.txt | awk '$2=="status"{print $1}')
for script in ${steplist}; do
     SLURM_ARRAY_TASK_ID=$statusIDNo bash $script $statusIDNo
done
