#!/bin/bash

# Marie-Madlen Pust

for i in `ls *.ndp.trm.s.mq30.mm.dup.calmd.bam.txt.rep.us.csv`

        do
                echo $i
                sbatch -c 16 combine_raspir_wochenende.sh $i
done