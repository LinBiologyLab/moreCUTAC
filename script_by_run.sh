runIDNo=$1
steplist=$(grep -v '#' steplist.txt | awk '$2=="run"{print $1}')
for script in ${steplist}; do
    SLURM_ARRAY_TASK_ID=$runIDNo bash $script $runIDNo
done
