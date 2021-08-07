
# README <br><br>
Combine the shotgun metagenomic sequencing alignment pipeline (Wochenende) with raspir <br>
For more information on Wochenende, see https://github.com/MHH-RCUG/Wochenende <br><br>


## Workflow 
Step 1. After running Wochenende (+ report), download the scripts
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
```

Step 5. Merge multiple files
```bash
sbatch multi_to_single.sh
```

Step 6. Rewrite id labels, add taxonomy information
```bash
python add_taxonomy.py wochenende.rep.us.raspir.merged.csv wochenende.rep.us.raspir.merged
```
Note: If the output file has missing taxonomy data, please contact me and provide the complete taxonomy details (domain, phylum, class, order, family, genus, and species) of the species you would like to have added. 
