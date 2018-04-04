#!/bin/sh
#SBATCH --job-name=lobSTR
#SBATCH --account=nn9244k
#SBATCH --time=12:0:0
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=16
#SBATCH --output=/PATH/TO/slurm_output/slurm-%j-%N.out

# I. Downloads reads from SRA stemming from the 1001 arabidopsis project (1135 accessions)
# II. Unpacks the reads to FASTA
# III. Runs lobSTR on the indexed TAIR10 genome
# IV. Deletes reads and FASTA files

# Loading modules
module load sratoolkit/2.7.0
module load lobstr/4.0.0b
module load samtools/1.3.1

# Setting variable paths
TMP=/PATH/TO/tmp/TRF_30
DATA=/PATH/TO/AccesionLists

# $1 comes from RunLobSTRMaster.sh, it is a list of six accessions
input=$1

# Moving to temporary folder
cd $TMP/data

# Looping over list and running lobSTR
while IFS= read -r var
do

echo "Downloading:" $var

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR194/$var/$var.sra

# Getting fasta files
fastq-dump.2.7.0 --split-files --fasta 0 $var.sra

# Running lobSTR (v. 3.0.2)

lobSTR --p1 "$var"_1.fasta --p2 "$var"_2.fasta \
--index-prefix /PATH/TO/GCA_000001735.1.TAIR10.12345_index/lobSTR_ \
-o $var \
-p 16 \
--rg-sample $var \
--rg-lib $var

# Sort and index using samtools:

samtools sort $var.aligned.bam -o $var.sorted.bam
samtools index $var.sorted.bam

# Delete RAW data to save space!
rm $var.sra
rm "$var"_1.fasta
rm "$var"_2.fasta

done < $DATA/$input