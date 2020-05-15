#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=03:59:00,mem=4gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

# 1. Extract the SNPs from the call set
java -Xmx4g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ${Var} \
-selectType SNP \
-o ./data/hard_filtered_variants/collective_samples_raw_snps_extracted.g.vcf \
2> ./analysis/hard_filtered_variants/collective_samples_raw_snps_extracted.log

if test -f ./data/hard_filtered_variants/collective_samples_raw_snps_extracted.g.vcf
then
  echo "Raw SNP Extraction Complete: ./data/hard_filtered_variants/collective_samples_raw_snps_extracted.g.vcf $(date +%F)$" >> \
  ./analysis/hard_filtered_variants/collective_samples_raw_snps_extracted.log
fi

print "fin"
