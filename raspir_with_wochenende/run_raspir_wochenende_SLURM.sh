#!/bin/bash

# run_raspir_wochenende_SLURM.sh
# last updated: 04 August 2021

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 30000

# set run x cpus
#SBATCH --cpus-per-task 8

# set name of job
#SBATCH --job-name=raspir_wochenende_run

# Add miniconda3 to PATH. TODO - autodetection
. /mnt/sfb900nfs/groups/tuemmler/mariep/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate raspir_env >> /dev/null

input_csv=$1
output_prefix=${input_csv%.ndp.trm.s.mm.dup.mq30.raspir.csv}

echo $input_csv
echo $output_prefix

srun python raspir_wochenende.py $input_csv $output_prefix
