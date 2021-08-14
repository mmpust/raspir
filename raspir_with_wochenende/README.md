# README <br><br>
Combine the shotgun metagenomic sequencing alignment pipeline (Wochenende) with raspir <br>
For more information on Wochenende, see https://github.com/MHH-RCUG/Wochenende <br><br>


## Workflow 
Step 1. After running Wochenende (+ report), download the scripts for running raspir and for combining Wochenende and raspir outputs
```bash
git clone https://github.com/mmpust/raspir.git

# copy unsorted Wochenende reporting files to the new directory
cp *.ndp.trm.s.mm.dup.mq30.calmd.bam.txt.rep.us.csv raspir/raspir_with_wochenende/   
# copy .BAM files for use with raspir 
cp *.ndp.trm.s.mm.dup.mq30.bam raspir/raspir_with_wochenende/   

# go to directory
cd raspir/raspir_with_wochenende/
```

Step 2. Prepare input files for raspir run from Wochenende files (.BAM):
```bash
sbatch prepare_file_for_raspir.sh   # adjust your conda environment at the top of the script
```

Step 3. Run raspir on .CSV files
```bash
sbatch runbatch_raspir_wochenende_SLURM.sh   # adjust your conda environment at the top of the script
```

Step 4. Combine raspir and wochenende output
```bash
sbatch runbatch_combine_raspir_wochenende_SLURM.sh   # adjust your conda environment at the top of the script

# Note: By default, the script will extract the normalised read counts (bacterial cell to human cell ratio). 
# If you want the raw read counts, open the script and change the column information:
nano combine_raspir_wochenende.sh
# go to line 36, change argument $8 to $3
# original file
awk -F',' '{print $1","$8}' ${input_file_wochenende_short}.temp_1.csv > ${input_file_wochenende_short}.temp_2.csv
# new file
awk -F',' '{print $1","$3}' ${input_file_wochenende_short}.temp_1.csv > ${input_file_wochenende_short}.temp_2.csv
# run file
sbatch runbatch_combine_raspir_wochenende_SLURM.sh 
```

Step 5. Merge multiple files
```bash
sbatch multi_to_single.sh
```

Step 6. Rewrite id labels, add taxonomy information
```bash
python add_taxonomy.py wochenende.rep.us.raspir.merged.csv wochenende.rep.us.raspir.merged
```

A table is generated (.CSV format). The assignment output has 7 taxonomy-specific columns. <br>

Domain | Phylum  | Class  | Order | Family | Genus | Species | Sample_1 | Sample_2 
---   | --- | --- | --- | --- | ---  | --- | --- | ---
Bacteria | Actinobacteria | Actinobacteria | Bifidobacteriales | Bifidobacteriaceae | Bifidobacterium | Bifidobacterium bifidum | 2.01 | 0.00
Bacteria | Proteobacteria | Gammaproteobacteria | Pasteurellales | Pasteurellaceae | Haemophilus | Haemophilus parainfluenzae | 40.1 | 52.3

<br>Note: If the output file has missing taxonomy data, please contact me and provide the complete taxonomy details (domain, phylum, class, order, family, genus, and species) of the species you would like to have added. 

