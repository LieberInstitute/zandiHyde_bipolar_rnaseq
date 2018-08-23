#!/bin/bash

## Usage:
# sh compute_mean.sh

mkdir -p logs

for region in amygdala sacc dlpfc
do
    for casestatus in "all" "case" "control"
    do
        
        SHORT="mean_${region}_${casestatus}"

# Construct shell file
echo "Creating script for ${SHORT}"
cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=72G,h_vmem=72G,h_fsize=200G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

module load ucsctools
module load wiggletools/default
Rscript compute_mean.R -r ${region} -c ${casestatus}

echo "**** Job ends ****"
date
EOF
        call="qsub .${SHORT}.sh"
        echo $call
        $call
    done
done
