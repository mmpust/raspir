# README
# Combine the output of the shotgun metagenomic sequencing alignment pipeline (Wochenende) with the raspir output
# For more information on Wochenende, see https://github.com/MHH-RCUG/Wochenende

# Step 1.
# run Wochenende and Wochenende report
# write the following new directories directly to the Wochenende output
```bash
mkdir run_raspir/ 
mkdir raspir_wochenende/

cp *.rep.us.* raspir_wochenende/ # copy unsorted Wochenende reporting files to the new directory
cp *.bam run_raspir/ # copy .BAM files for use with raspir 
```

# Step 2.
# download the scripts
```bash
git clone .
cp raspir_with_wochenende/*.sh .
cp raspir_with_wochenende/*.py .
cp raspir_wochenende.py run_raspir_Wochenende_SLURM.sh runbatch_raspir_wochenende_SLURM.sh run_raspir/
cp combine_raspir_wochenende.sh multi_to_single_file.sh runbatch_combine_raspir_wochenende_SLURMs.sh run_taxonomy.py taxonomy_file raspir_wochenende/
```

# prepare input files for raspir run from Wochenende bam files:
```bash
cd run_raspir/
sbatch prepare_file_for_raspir.sh # insert your conda environment at the top of the script
```

# run raspir on .CSV files
```bash
sbatch runbatch_raspir_wochenende_SLURM.sh
```

# go back and copy raspir output into raspir_wochenende/ directory
```bash
cd ..
cp run_raspir/*.csv raspir_wochenende/
cd raspir_wochenende/
```

# combine raspir and wochenende output
```bash
sbatch runbatch_combine_raspir_wochenende_SLURM.sh
```

# create one single dataframe from multiple files
```bash
sbatch multi_to_single_SLURM.sh
```

# rewrite id labels, add taxonomy information
```bash
python add_taxonomy.py wochenende.rep.us.raspir.merged.csv
```
