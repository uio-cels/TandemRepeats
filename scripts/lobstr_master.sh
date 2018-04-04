#!/bin/sh

# What it does:
# I. Splits the list of 1135 reads into lists of three reads per file
# II. Starts lobSTR_Slave.slurm that runs lobSTR on the reads

# Setting variables

SRALIST=1135_SraAccList.txt # List of with accession IDs
DATA=/PATH/TO/AccesionLists

# Loading modules
module load sratoolkit/2.7.0

# Splitting list into lists of six

cd $DATA

split -l 6 $SRALIST SRA

# Starting the jobs doing the work

for LIST in $(ls | grep "SRA");

do
echo $LIST
sbatch /PATH/TO/lobstr_slave.slurm $LIST

done