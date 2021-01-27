#!/bin/bash

## Marie-Madlen Pust
## Last updated: 27 January 2021
## pust.marie-madlen@mh-hannover.de

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
#SBATCH --cpus-per-task 12

# set name of job
#SBATCH --job-name=raspir_prepare

# Add miniconda3 to PATH
. /miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate raspir_env >> /dev/null

for items in *.sam
	do
		echo $items
		fname=$(echo ${items} | sed 's/.sam//')
		echo $fname

    		# Remove reads with low mapping quality
    		samtools view -hM -q 20 $items > ${items%.sam}.mq20.sam

    		# Convert file from SAM to BAM format
    		samtools view -h -b -S ${items%.sam}.mq20.sam  > ${fname}.bam

    		# Discard unmapped sequences
    		samtools view -b -F 4 $items > ${fname}_1.bam

    		# Sort bam file
    		samtools sort ${fname}_1.bam -o ${fname}.sorted.bam

    		# Generate index
    		samtools index ${fname}.sorted.bam

    		# Obtain coverage information
    		samtools depth ${fname}.sorted.bam > ${fname}.sorted.raspir1.csv

    		# Add genome size
    		sed 's/LN://g' ${items} > ${fname}.genomeSize_1.csv
    		sed -i 's/SN://g' ${fname}.genomeSize_1.csv
    		cut -f2- ${fname}.genomeSize_1.csv > ${fname}.genomeSize.csv

    		# Add column with genome size
    		awk -v FS="\t" -v OFS="\t" 'FNR==NR{a[$1]=$2;next;} {if(a[$1]) {print a[$1], $0} else {print "NA",$0}}' \
		${fname}.genomeSize.csv ${fname}.sorted.raspir1.csv > ${fname}.sorted.raspir.csv

		# Convert into a comma-separated file
		sed -i 's/\t/,/g' ${fname}.sorted.raspir.csv
		# Add header
		sed -i '1iGenomeLength,Organism,Position,Depth\' ${fname}.sorted.raspir.csv

		# Remove intermediate files
		rm ${fname}.sorted.raspir1.csv ${fname}.genomeSize_1.csv ${fname}.genomeSize.csv
	done
