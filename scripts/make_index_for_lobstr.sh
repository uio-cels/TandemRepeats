#!/bin/sh

GENOME=/PATH/TO/GCA_000001735.1.TAIR10.12345.fasta

module load bedtools/2.25.0
module load python2/2.7.10
module load lobstr/4.0.0b

python /PATH/TO/convert_trf_bed_lobstr.py -d /PATH/TO/GCA_000001735.1.TAIR10.12345.fasta.2.7.7.80.10.30.3.dat > GCA_000001735.1.TAIR10.12345.fasta_trf.bed

mkdir -p GCA_000001735.1.TAIR10.12345_index
lobstr_index.py --str GCA_000001735.1.TAIR10.12345.fasta_trf.bed --ref $GENOME --out GCA_000001735.1.TAIR10.12345_index 1> lobstr_index.out 2> lobstr_index.err

GetSTRInfo.py GCA_000001735.1.TAIR10.12345.fasta_trf.bed $GENOME > GCA_000001735.1.TAIR10.12345_index/GCA_000001735.1.TAIR10.12345_strinfo.tab