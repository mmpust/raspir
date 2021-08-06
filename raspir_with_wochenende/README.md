
# README <br><br>
Combine the output of the shotgun metagenomic sequencing alignment pipeline (Wochenende) with raspir <br>
For more information on Wochenende, see https://github.com/MHH-RCUG/Wochenende <br><br>


## Workflow 
Step 1. After running Wochenende (+ report), download the scripts
```bash
git clone .
cp *.rep.us.* raspir_with_wochenende/ # copy unsorted Wochenende reporting files to the new directory
cp *.bam raspir_with_wochenende/ # copy .BAM files for use with raspir 
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
```

Step 5. Merge multiple files
```bash
sbatch multi_to_single.sh
```

Step 6. Rewrite id labels, add taxonomy information
```bash
python add_taxonomy.py wochenende.rep.us.raspir.merged.csv
```
