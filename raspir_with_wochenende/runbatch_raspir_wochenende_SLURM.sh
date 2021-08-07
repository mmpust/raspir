#!/bin/bash

# Marie Pust
# last updated: 04 August 2021

for i in `ls *.ndp.trm.s.mm.dup.mq30.raspir.csv`
	do
		echo $i
		sbatch -c 16 run_raspir_wochenende_SLURM.sh $i
done
