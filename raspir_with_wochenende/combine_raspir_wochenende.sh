#!/bin/bash
# Author: Marie-Madlen Pust
# Written: 04 August 2021

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 30000

# set run x cpus
#SBATCH --cpus-per-task 8

# set name of job
#SBATCH --job-name=combine_raspir_wochenende

# Add miniconda3 to PATH. TODO - autodetection
. /mnt/sfb900nfs/groups/tuemmler/mariep/miniconda3/etc/profile.d/conda.sh

# merge wochenende and raspir output
# unsorted file (.CSV format), after running basic_reporting.py
input_file_wochenende=$1

# do not edit!
input_file_wochenende_short=${input_file_wochenende%.ndp.%}
input_file_raspir=${input_file_wochenende_short}.final_stats.csv

echo $input_file_wochenende
echo $input_file_wochenende_short
echo $input_file_raspir

# remove human reads from wochenende file
sed '/^1_1_1_/d' $input_file_wochenende > ${input_file_wochenende_short}.temp_1.csv

# returns only taxon id and normalised read count (normalised to human reads) from wochenende output
awk -F',' '{print $1","$8}' ${input_file_wochenende_short}.temp_1.csv > ${input_file_wochenende_short}.temp_2.csv

# returns only taxon id and uniform/non-uniform column from raspir output
awk -F',' '{print $1","$6}' $input_file_raspir > ${input_file_raspir}.temp_1.csv

# if taxon id (raspir output) matches with id (wochenende output)
join -t, <(sort ${input_file_wochenende_short}.temp_2.csv) <(sort ${input_file_raspir}.temp_1.csv) -a 1 -o auto -e 'non_uniform' > ${input_file_wochenende_short}.temp_2.raspir.csv

# if prediction is "non-uniform", set count value to 0
awk -F',' '{ if($3=="uniform") $4=$2; else $4=0; print $0; } ' ${input_file_wochenende_short}.temp_2.raspir.csv > ${input_file_wochenende_short}.temp_3.raspir.csv

# remove old columns
awk '{$2=$3=""; print $0}' ${input_file_wochenende_short}.temp_3.raspir.csv > ${input_file_wochenende_short}.temp_4.raspir.csv

# sort taxon id column alphabetically
sort -t',' -k1 ${input_file_wochenende_short}.temp_4.raspir.csv > ${input_file_wochenende_short}.temp_5.raspir.csv

# add header
sed -i 1i"organism,bcphc_normalised" ${input_file_wochenende_short}.temp_5.raspir.csv

# replace tab with comma
sed -e 's/   /,/g' ${input_file_wochenende_short}.temp_5.raspir.csv > ${input_file_wochenende_short}.wochenende.rep.us.raspir.csv
