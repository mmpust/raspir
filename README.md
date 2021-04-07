![Alt text](raspir_image.jpg?raw=true)
# Background
In shotgun metagenomic sequencing experiments, the total DNA is extracted from complex samples. The DNA of rare species is scarce, and poor coverage of these genomes is often observed. Scientists commonly define thresholds to exclude 1-10% of the least abundant species in a given sample from further analyses. On the one hand, this filtering step allows for robust investigations of core communities. On the other hand, valuable information on the community structure is lost. The rare biosphere harbours more species than the core microbial community and hence provides the environment of interest with high functional flexibility.

However, if only a few short DNA reads are detected that are unique to species A, there are at least three explanations: <br>
a) Sample contamination; <br>
b) Rare species A was present in the environment of interest, so it is a true positive species. In this case, the reads are expected to spread across the entire reference genome in a fairly uniform manner due to random DNA sequencing; or <br>
c) Rare species A was absent (false positive) but rare species B was present, which acquired genes of species A during past events. In this case, the reads are expected to cluster at specific sides of the reference genome of species A. <br>

The raspir tool calculates a position-domain signal based on the distances of reads aligning to a circular reference genome and converts the information into a frequency signal through Discrete Fourier Transforms (DFT). In addition, a reference frequency signal is constructed with the same number of reads, but with an ideal uniform distribution of reads across the reference genome. Both frequency signals are compared using Pearson's correlation measures. 

# Implementation
Raspir is implemented in Python 3.7. Using this tool, a real-world dataset (5 GB) containing information on hundreds of species can be processed on a single node server. The input data must be structured in the following manner: genome length of the corresponding reference genome, organism, read position, read depth. 

See the following section for further information on how to convert .FASTQ files into the .CSV input files to sucessfully execute raspir.

# Get started
## Set up the environment
### Install conda packages 

```bash
conda create --name raspir_env
conda activate raspir_env

# Install trimmomatic [1] 
conda install -c bioconda trimmomatic
# Install samtools [2] 
conda install -c bioconda samtools

# Install your alignment tool of choice
# Burrows-Wheeler aligner [3]
conda install -c bioconda bwa
# Bowtie2 [4]
conda install -c bioconda/label/cf201901 bowtie2

# Install python packages for raspir
conda install pandas
conda install -c conda-forge statsmodels
conda install -c conda-forge matplotlib
```


### Create working directory

```bash
mkdir raspir/
cd raspir/

YOURPATH=$PWD
echo "$YOURPATH"

mkdir reference_database/
mkdir run_raspir/
```

### Set-up the reference database <br>
The reference database has to be downloaded and indexed only once. 
The database contains only complete bacterial genomes with one strain per species.
Note: You may also use a customised reference database. It is however strongly recommended to avoid draft or low-quality reference genomes and use complete sequences of circular microorganisms only.

```bash
# Load database into your working directory
cd reference_database/
# Get database from
https://sync.academiccloud.de/index.php/s/wSFAaUpfkjY4fZG/download
# Unzip the reference database 
gunzip complete_bacterialRefSeqs_201910_3.fasta

# Generate an index of the reference fasta depending on the alignment tool of your choice
samtools faidx complete_bacterialRefSeqs_201910_3.fasta
bwa index complete_bacterialRefSeqs_201910_3.fasta
bowtie2-build complete_bacterialRefSeqs_201910_3.fasta complete_bacterialRefSeqs_201910_3
cd ..
```

## Data cleaning, alignment & sorting
### Trimmomatic 
Load your .FASTQ files into the working directory and run Trimmomatic for quality trimming and adapter clipping.

```bash
cd run_raspir/

# Paired-end data
trimmomatic PE \
  R1.fastq R2.fastq \
  R1.trim.fastq R1.un.trim.fastq \
  R2.trim.fastq R2.un.trim.fastq \
  SLIDINGWINDOW:4:20 MINLEN:25

# Single-end reads
trimmomatic SE \
  R.fastq \
  R.trim.fastq \
  SLIDINGWINDOW:4:20 MINLEN:25
```

### Alignment 
```bash
# Burrows-Wheeler-Aligner (bwa)
# see http://bio-bwa.sourceforge.net/
# For paired-end reads
bwa mem $YOURPATH/reference_database/complete_bacterialRefSeqs_201910_3.fasta \
  R1.trim.fastq R2.trim.fastq > R.trim.bwa.sam
# For single-end reads
bwa mem $YOURPATH/reference_database/complete_bacterialRefSeqs_201910_3.fasta \
  R.trim.fastq > R.trim.bwa.sam

# Bowtie2
# see http://bowtie-bio.sourceforge.net/
# For paired-end reads
bowtie2 -x $YOURPATH/reference_database/complete_bacterialRefSeqs_201910_3 \
  -1 R1.trim.fastq -2 R2.trim.fastq -S R.trim.bowtie2.sam
# For single-end reads
bowtie2 -x $YOURPATH/reference_database/complete_bacterialRefSeqs_201910_3 \
  -U R.trim.fastq -S R.trim.bowtie2.sam
```

### Sorting, indexing & final clean-up

```bash
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
```

# Run raspir
Download the python script into the run_raspir/ folder.
If you have multiple input .CSV files in your directory, the script will automatically create a list of input files and iterate over the list.

```python
# Linux
python3 raspir_v1_0_unix.py

# Windows
python3 raspir_v1_0_win.py
```

# Output
A table is generated (.CSV format). The assignment output has 6 columns.

| Species | r_value  | p_value  | stError | euclidean | distribution |
| :---:   | :-: | :-: | :-: | :-: | :-: | 
| Pseudomonas aeruginosa | 0.99 | 0.0 | 0.00019 | 0.01 | uniform |
| Streptococcus salivarius | 0.97 | 0.0 | 0.00016 | 0.002 | uniform |
| Rothia mucilaginosa | 0.99 | 0.0 | 0.000002 | 0.0001 | uniform | 


# Contributors
@mmpust author
@colindaven updates
@nick-youngblut updates 

# References
[1] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170. <br>
[2] Li H., Handsaker B., Wysoker A. et al. (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. <br>
[3] Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. <br>
[4] Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359. <br>
