# README



# Step 1: Run Wochenende (https://github.com/MHH-RCUG/Wochenende)
mkdir run_raspir/ 
mkdir raspir_wochenende/

cp *.us.rp.csv raspir_wochenende/ # copy unsorted Wochenende reporting files to the new directory
cp *.bam run_raspir/ # copy .BAM files to the new directory 

# prepare input files for raspir run from Wochenende bam files:
cd run_raspir/

# run raspir on .CSV files

# go back and copy raspir output into raspir_wochenende/ directory
cd ..
cp run_raspir/*.csv raspir_wochenende/
cd raspir_wochenende/



/ngsssd1/tuem_mp/projects_MariePust_2018/longitudinial_cf_analysis/run_centrifuge

/mnt/sfb900nfs/groups/tuemmler/mariep/programs/wochenende
2021_08_meta_fungi_human_masked
/mnt/sfb900nfs/groups/tuemmler/mariep/programs/wochenende
