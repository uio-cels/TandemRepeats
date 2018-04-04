#!/bin/sh

cd /PATH/TO/data # Containing all the BAM files

ls | grep "sorted.bam" |grep -v ".bai" > AllGenesAllAccessions.bamlist

#	Loading bamtools
module load bamtools/2.3.0

#	Running bamtools

LINES=$(wc -l AllGenesAllAccessions.bamlist | awk '{print $1}')
HALF=$(echo $(($LINES / 2)))

# Splitting file into two, as bamtools cannot take more than 1024 files
split -l $HALF AllGenesAllAccessions.bamlist AllGenesAllAccessions.bamlist_

# 	Running first half
bamtools merge -list AllGenesAllAccessions.bamlist_aa -out AllGenesAllAccessions_aa.merged.bam
bamtools index -in AllGenesAllAccessions_aa.merged.bam

#	Running second half
bamtools merge -list AllGenesAllAccessions.bamlist_ab -out AllGenesAllAccessions_ab.merged.bam
bamtools index -in AllGenesAllAccessions_ab.merged.bam

#	Merging the merged 
bamtools merge -in AllGenesAllAccessions_aa.merged.bam \
-in AllGenesAllAccessions_ab.merged.bam \
-out AllGenesAllAccessions.merged.bam

bamtools index -in AllGenesAllAccessions.merged.bam


#	Loading lobSTR, contains allelotype
module load lobstr/4.0.0b

# Running allelotype
allelotype \
  --command classify \
  --bam AllGenesAllAccessions.merged.bam \
  --noise_model /PATH/TO/lobSTR/models/illumina_v2.0.3 \
  --out /PATH/TO/allelotype/AllGenesAllAccessions.merged.bam \
  --strinfo /PATH/TO/GCA_000001735.1.TAIR10.12345_index/GCA_000001735.1.TAIR10.12345_strinfo.tab \
  --index-prefix /PATH/TO/GCA_000001735.1.TAIR10.12345_index/lobSTR_ \
  > allelotype.out 2> allelotype.err

# Annotate

module load tabix/0.2.6
module load bcftools/1.3
module load vcftools/0.1.11

# First sorting
vcf-sort AllGenesAllAccessions.merged.bam.vcf > AllGenesAllAccessions.merged.bam.sorted.vcf

# NB! REMOVE THE ^M FROM THE sv_gene.data file!

## sed -e "s/^M//" sv_gene.data > sv_gene.data_nom, but need to write  ^M in the command line!

# Then making a bed file in the correct format
awk '{print "chr" $6 "\t" $3 "\t" $4 "\t" $2}' sv_gene.data_nom > sv_gene.data.bed


# Zipping and tabbing bed
bgzip -c sv_gene.data.bed > sv_gene.data.bed.gz
tabix -p bed sv_gene.data.bed.gz

# Zipping and tabbing vcf
bgzip -c AllGenesAllAccessions.merged.bam.sorted.vcf > AllGenesAllAccessions.merged.bam.sorted.vcf.gz
tabix -p vcf AllGenesAllAccessions.merged.bam.sorted.vcf.gz

echo "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">
##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr5>" > header_file.txt

# Annotating
bcftools annotate \
-a sv_gene.data.bed.gz \
-c CHROM,FROM,TO,GENE \
-h header_file.txt \
-o AllGenesAllAccessions.merged.bam.sorted.annotated.vcf \
AllGenesAllAccessions.merged.bam.sorted.vcf