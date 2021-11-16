#!/bin/bash

# Marie Pust
# last updated: 04 August 2021

for i in `ls *.raspir.csv`
	do
		echo $i
		sbatch -c 16 run_raspir_wochenende_SLURM.sh $i
done
