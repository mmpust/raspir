
# README <br><br>
Combine the output of the shotgun metagenomic sequencing alignment pipeline (Wochenende) with raspir <br>
For more information on Wochenende, see https://github.com/MHH-RCUG/Wochenende <br><br>


## Workflow 
Step 1. After running Wochenende (+ report), create the following two directories in the same folder
```bash
mkdir run_raspir/ 
mkdir raspir_wochenende/

cp *.rep.us.* raspir_wochenende/ # copy unsorted Wochenende reporting files to the new directory
cp *.bam run_raspir/ # copy .BAM files for use with raspir 
```

Step 2. Download the scripts
```bash
git clone .
cp raspir_with_wochenende/*.sh .
cp raspir_with_wochenende/*.py .
cp raspir_wochenende.py run_raspir_Wochenende_SLURM.sh runbatch_raspir_wochenende_SLURM.sh run_raspir/
cp combine_raspir_wochenende.sh multi_to_single_file.sh runbatch_combine_raspir_wochenende_SLURMs.sh run_taxonomy.py taxonomy_file raspir_wochenende/
```

Step 3. Prepare input files for raspir run from Wochenende files (.BAM):
```bash
cd run_raspir/
sbatch prepare_file_for_raspir.sh # adjust your conda environment at the top of the script
```

Step 4. Run raspir on .CSV files
```bash
sbatch runbatch_raspir_wochenende_SLURM.sh
```

Step 5. Combine raspir and wochenende output
```bash
cd ..
cp run_raspir/*.csv raspir_wochenende/
cd raspir_wochenende/

sbatch runbatch_combine_raspir_wochenende_SLURM.sh
```

Step 6. Create a single dataframe from multiple files
```bash
sbatch multi_to_single_SLURM.sh
```

Step 7. Rewrite id labels, add taxonomy information
```bash
python add_taxonomy.py wochenende.rep.us.raspir.merged.csv
```
