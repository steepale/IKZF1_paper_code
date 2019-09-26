#===============================================================================
#
#         FILE: /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/GATK_Variant_Analysis_Germline_DNA/scripts/GATK_Variant_Analysis_Germline_DNA_Main_Documentation.sh
#
#        USAGE: ./scripts/GATK_Variant_Analysis_Germline_DNA_Main_Documentation.sh
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for Germline variant
#                calls (SNPs and Indels) on preprocessed BAM files
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2016.05.18
#     REVISION:  
#===============================================================================

# Permanent PROJECT DIRECTORYs (MSU Cluster)
# /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels
# /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/GATK_Variant_Analysis_Germline_DNA
# Scratch PROJECT DIRECTORY (MSU Cluster)
# /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

if [ "$#" -ne 2 -o (! -r "$1" -a ! -r "$2") ]
then
	echo "Inappropriate Input File Amount or Input Files are not Readable"
	exit 1
fi

# Calling germline variants with GATK
# https://www.broadinstitute.org/gatk/guide/article?id=2803

# The majority of analysis will be performed in scratch and the neccessary files will be moved to the permanent project directory
cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# make appopriate directories
mkdir -p ./{data,scripts,analysis}

mkdir ./data/ref
mkdir ./data/bams
mkdir ./data/raw_snps_indels
mkdir ./analysis/raw_snps_indels
mkdir ./data/hard_filtered_variants
mkdir ./analysis/hard_filtered_variants
mkdir ./analysis/hard_filtered_variants/R_plots
mkdir ./data/truth_files
mkdir ./analysis/liftover
mkdir ./data/truth_files/divide_ma
mkdir ./scripts/compare_snps
mkdir ./data/vsqr
mkdir ./scripts/picard
mkdir ./data/truth_files/divide_ma/final_vcf
mkdir ./scripts/snpEff
mkdir ./data/germline_snps
mkdir ./analysis/snpEff
mkdir ./analysis/vsqr
mkdir ./data/germline_snps/indv_samples
mkdir ./scripts/indv_samples
mkdir ./data/arrays
mkdir ./data/germline_snps/indv_samples/tumors
mkdir ./data/germline_snps/indv_samples/normals
mkdir ./data/germline_snps/indv_samples/parents
mkdir ./data/germline_snps/indv_samples/private_snps_divided
mkdir ./data/haplosaurus
mkdir ./data/transmission_phased
mkdir ./data/GenVisR
mkdir ./data/proteinseqs

# Create symbolic links of all the BAM files to the data directory.
ln -s "/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/final_bam/"*"_Bwa_RG_dedupped_realigned.bam" \
./data/bams/

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# Call Germline raw variants (SNPs and Indels) on DNA Sequencing samples
find ./data/bams -name "*-0_S*_Bwa_RG_dedupped_realigned.bam" |\
xargs -i echo 'qsub ./scripts/GATK_haplotype_caller_raw_snps_indels_GVCF.sh -v Var='{} |sh

# ./scripts/GATK_haplotype_caller_raw_snps_indels_GVCF.sh
################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=06:00:00:00,mem=120gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

sample_name=$(basename ${Var} "_Bwa_RG_dedupped_realigned.bam")

java -Xmx120g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ./data/ref/galgal5.fa \
-I ${Var} \
-ERC GVCF \
-o ./data/raw_snps_indels/${sample_name}_raw_snps_indels.g.vcf \
2> ./analysis/raw_snps_indels/${sample_name}_raw_snps_indels.log

if test -f ./data/raw_snps_indels/${sample_name}_raw_snps_indels.g.vcf
then
	echo "Raw Germline Variants Call Complete: ./data/raw_snps_indels/${sample_name}_raw_snps_indels.g.vcf $(date +%F)$" >> \
  ./analysis/raw_snps_indels/${sample_name}_raw_snps_indels.log
fi
##################################

# Perform Joint Genotyping
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

# Make sure th files are in uncompressed format if need be
# find ./data/raw_snps_indels -name "*.g.vcf.gz" |
# xargs -i sh -c "gunzip {}"

# Compress the files with bgzip (this process will take very long and you may need to edit the script)
module load tabix/0.2.6

find ./data/raw_snps_indels/ -name "[0-9]*-0*S*.g.vcf" | \
xargs -i sh -c "bgzip {}"

# Index the vcf files
find ./data/raw_snps_indels/ -name "[0-9]*-0*S*.g.vcf.gz" | \
xargs -i sh -c "tabix -p vcf {}"


# Write the joint genotyping script before it can be submitted
cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels
cat <<Header_input > ./scripts/joint_genotype_germline_gvcfs.sh
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=72:00:00,mem=60gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

java -Xmx60g -cp \$GATK -jar \$GATK/GenomeAnalysisTK.jar \\
-T GenotypeGVCFs \\
-R ./data/ref/galgal5.fa \\
Header_input

find ./data/raw_snps_indels/ -name "[0-9]*-0*S*.g.vcf.gz" | \
xargs -i printf '%s\n' '-V {} \' >> ./scripts/joint_genotype_germline_gvcfs.sh

cat <<Footer_input >> ./scripts/joint_genotype_germline_gvcfs.sh
-o ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz \\
2> ./analysis/raw_snps_indels/germline_raw_snps_indels_genotyped.log

if test -f ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz
then
  echo "Joint Genotyping Complete: ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz \$(date +%F)\$" >> \\
  ./analysis/raw_snps_indels/germline_raw_snps_indels_genotyped.log
fi

qstat -f \${PBS_JOBID}
Footer_input

# # Perform Joint Genotyping on multiple samples
qsub ./scripts/joint_genotype_germline_gvcfs.sh

# Perform Joint genotyping on a single sample:
# Line 6
find `pwd` -name "002683_Line-6_raw_snps_indels.g.vcf" |\
xargs -i echo 'qsub ./scripts/joint_genotype_single_sample.sh -v Var='{} |sh

# Line 7
find `pwd` -name "002684_Line-7_raw_snps_indels.g.vcf" |\
xargs -i echo 'qsub ./scripts/joint_genotype_single_sample.sh -v Var='{} |sh

# ./scripts/joint_genotype_single_sample.sh
################################
# #!/bin/bash -login
# #PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=60gb
# #PBS -j oe
# set -e
# set -u
# set -o pipefail
# 
# module load GATK/3.5.0
# 
# cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels
# 
# sample_name=$(basename ${Var} "_raw_snps_indels.g.vcf")
# 
# java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
# -T GenotypeGVCFs \
# -R ./data/ref/galgal5.fa \
# -V ${Var} \
# -o ./data/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.g.vcf \
# 2> ./analysis/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.log
# 
# if test -f ./data/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.g.vcf
# then
#   echo "Joint Genotyping Complete: ./data/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.g.vcf $(date +%F)$" >> \
#   ./analysis/raw_snps_indels/${sample_name}_raw_snps_indels_genotyped.log
# fi
# 
# qstat -f ${PBS_JOBID}
##################################

###########################################################################################################

# Apply hard filters to the call set 
# https://www.broadinstitute.org/gatk/guide/article?id=2806

# We apply hard filters very stringently
# The strignent calls will account for the confident Germline variants already present 

# https://www.broadinstitute.org/gatk/guide/article?id=2806
# 1. Extract the SNPs from the call set
# 2. Determine parameters for filtering SNPs
# 3. Apply the filter to the SNP call set
# 4. Extract the Indels from the call set
# 5. Determine parameters for filtering indels
# 6. Apply the filter to the Indel call set

# 1. Extract the SNPs from the call set

# line 6x7 F1
find ./data/raw_snps_indels -name "germline_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/extract_snps_from_genotyper.sh -v Var='{} |sh

# Compress and index the parental genotyped vcf files of SNPs and indels
module load tabix/0.2.6

# Compress the vcf file
find ./data/raw_snps_indels/ -name "002683_Line-6_raw_snps_indels_genotyped.g.vcf" | \
xargs -i sh -c "bgzip {}"

find ./data/raw_snps_indels/ -name "002684_Line-7_raw_snps_indels_genotyped.g.vcf" | \
xargs -i sh -c "bgzip {}"

# Index the vcf files
find ./data/raw_snps_indels/ -name "002683_Line-6_raw_snps_indels_genotyped.g.vcf.gz" | \
xargs -i sh -c "tabix -p vcf {}"

find ./data/raw_snps_indels/ -name "002684_Line-7_raw_snps_indels_genotyped.g.vcf.gz" | \
xargs -i sh -c "tabix -p vcf {}"


# line 6
find ./data/raw_snps_indels/ -name "002683_Line-6_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/extract_snps_from_genotyper.sh -v Var='{} |sh

# line 7
find ./data/raw_snps_indels/ -name "002684_Line-7_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/extract_snps_from_genotyper.sh -v Var='{} |sh

# ./scripts/extract_snps_from_genotyper.sh
################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=03:59:00,mem=4gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

sample_name=$(basename ${Var} "_raw_snps_indels_genotyped.g.vcf.gz")

# 1. Extract the SNPs from the call set
java -Xmx4g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ${Var} \
-selectType SNP \
-o ./data/raw_snps_indels/${sample_name}_raw_snps_extracted.g.vcf.gz \
2> ./analysis/raw_snps_indels/${sample_name}_raw_snps_extracted.log

if test -f ./data/raw_snps_indels/${sample_name}_raw_snps_extracted.g.vcf.gz
then
  echo "Raw SNP Extraction Complete: ./data/raw_snps_indels/${sample_name}_raw_snps_extracted.g.vcf.gz $(date +%F)" >> \
  ./analysis/raw_snps_indels/${sample_name}_raw_snps_extracted.log
fi

print "fin"
#################################

# Count all the SNPs per sample to make sure sample calls worked
cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

find ./data/raw_snps_indels/ -name "[0-9]*-0*S*.g.vcf.gz" | \
xargs -i basename {} | \
sed 's/.g.vcf.gz//' | \
sort | uniq | \
xargs -i echo 'qsub ./scripts/test_count_variants.sh -v Var='{} |sh

# ./scripts/test_count_variants.sh
##############################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=00:30:00,mem=8gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# Select individual samples and filter out filtered calls and all non-variants
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/raw_snps_indels/germline_raw_snps_extracted.g.vcf.gz \
-o ./data/raw_snps_indels/germline_raw_snps_${Var}.vcf.gz \
-sn ${Var} \
-env \
-ef
##############################

# Perform variant counts with bcftools
module load bcftools/1.2

find ./data/raw_snps_indels/ -name "*-0*S*.vcf.gz" | \
xargs -i basename {} | \
sed 's/.vcf.gz//' | \
sort | uniq | \
xargs -i bcftools stats ./data/raw_snps_indels/germline_raw_snps_{}.vcf.gz >> ./data/raw_snps_indels/counts.txt

# Print the SNP counts into a file
(echo "SNP Counts Before Hard Filters" > ./data/raw_snps_indels/SNP_counts_pre_filter.txt; \
grep -e "number of SNPs:" -e "./data/raw_snps_indels/germline_raw_snps" ./data/raw_snps_indels/counts.txt >> \
./data/raw_snps_indels/SNP_counts_pre_filter.txt)

# Uncompress for python scripts
bgzip -d -c ./data/raw_snps_indels/germline_raw_snps_indels.g.vcf.gz > ./data/raw_snps_indels/germline_raw_snps_indels.g.vcf

# 2. Determine parameters for filtering SNPs
# Extract QD Stats for R
python2.7 ./scripts/Collect_QD.py \
./data/raw_snps_indels/germline_raw_snps_indels.g.vcf \
./data/hard_filtered_variants/germline_raw_snps_extracted_QD.txt

# Extract FS Stats for R
python2.7 ./scripts/Collect_FS.py \
./data/raw_snps_indels/germline_raw_snps_indels.g.vcf \
./data/hard_filtered_variants/germline_raw_snps_extracted_FS.txt

# Extract MQ Stats for R
python2.7 ./scripts/Collect_MQ.py \
./data/raw_snps_indels/germline_raw_snps_indels.g.vcf \
./data/hard_filtered_variants/germline_raw_snps_extracted_MQ.txt

# Extract MQRankSum Stats for R
python2.7 ./scripts/Collect_MQRankSum.py \
./data/raw_snps_indels/germline_raw_snps_indels.g.vcf \
./data/hard_filtered_variants/germline_raw_snps_extracted_MQRankSum.txt

# Extract ReadPosRankSum Stats for R
python2.7 ./scripts/Collect_ReadPosRankSum.py \
./data/raw_snps_indels/germline_raw_snps_indels.g.vcf \
./data/hard_filtered_variants/germline_raw_snps_extracted_ReadPosRankSum.txt


module load R/3.2.0

Rscript --vanilla ./scripts/hard_filt_plots.R

# ./scripts/hard_filt_plots.R
##########################################
# library("ggplot2")
# 
# # In R, Generate the QD Plot
# samples_QD = read.table("./data/hard_filtered_variants/germline_raw_snps_extracted_QD.txt")
# qplot(V1, data = samples_QD, geom = "density") + xlab("QD") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_QD_density_snps.png")
# 
# # In R, Generate the FS Plot (x-axis is logged)
# samples_FS = read.table("./data/hard_filtered_variants/germline_raw_snps_extracted_FS.txt")
# qplot(V1, data = samples_FS, geom = "density", log="x") + xlab("FS") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_FS_density_snps.png")
# # Cutoff right after right peak
# 
# # In R, Generate the MQ Plot
# samples_MQ = read.table("./data/hard_filtered_variants/germline_raw_snps_extracted_MQ.txt")
# qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_MQ_density_snps.png")
# 
# # Zoom in on the MQ peak at 40
# qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(39.0, 41.0)
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_MQ_density_40_snps.png")
# 
# # Zoom in on the MQ peak at 60
# qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(59.0, 61.0)
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_MQ_density_60_snps.png")
# # Discard anything that is not 60
# 
# # In R, Generate the MQRankSum Plot
# samples_MQRankSum = read.table("./data/hard_filtered_variants/germline_raw_snps_extracted_MQRankSum.txt")
# qplot(V1, data = samples_MQRankSum, geom = "density") + xlab("MQRankSum") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_MQRankSum_density_snps.png")
# # Remove anything less than -2
# 
# # In R, Generate the ReadPosRankSum Plot
# samples_ReadPosRankSum = read.table("./data/hard_filtered_variants/germline_raw_snps_extracted_ReadPosRankSum.txt")
# qplot(V1, data = samples_ReadPosRankSum, geom = "density") + xlab("ReadPosRankSum") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/germline_ReadPosRankSum_density_snps.png")
##########################################


# 3. Apply the filter to the SNP call set
# https://www.broadinstitute.org/gatk/guide/article?id=6925
# https://www.broadinstitute.org/gatk/guide/article?id=2806

# Line 6x7 F1
find ./data/raw_snps_indels -name "germline_raw_snps_extracted.g.vcf.gz" | \
xargs -i echo 'qsub ./scripts/hard_filter_SNPs.sh -v Var='{} |sh

# line 6
find ./data/raw_snps_indels -name "002683_Line-6_raw_snps_extracted.g.vcf.gz" | \
xargs -i echo 'qsub ./scripts/hard_filter_SNPs.sh -v Var='{} |sh

# line 7
find ./data/raw_snps_indels -name "002684_Line-7_raw_snps_extracted.g.vcf.gz" | \
xargs -i echo 'qsub ./scripts/hard_filter_SNPs.sh -v Var='{} |sh


# ./scripts/hard_filter_SNPs.sh
################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=03:59:00,mem=8gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

sample_name=$(basename ${Var} "_raw_snps_extracted.g.vcf.gz")

# Perform stringent hard filtering
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/ref/galgal5.fa \
-V ${Var} \
--filterExpression "MQ < 60.0 || MQ > 60.0 || QD < 16.0 || FS > 0.778 || MQRankSum < 0.5" \
--filterName "SNP_HARD_FILTER" \
-o ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps_stringent.g.vcf.gz \
2> ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_snps_stringent.log

if test -f ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps_stringent.g.vcf.gz
then
  echo "Stringent Hard Filter Complete: /data/hard_filtered_variants/${sample_name}_hard_filtered_snps_stringent.g.vcf.gz $(date +%F)$" >> \
  ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_snps_stringent.log
fi

# Perform lenient hard filtering
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/ref/galgal5.fa \
-V ${Var} \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "SNP_HARD_FILTER" \
-o ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps_lenient.g.vcf.gz \
2> ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_snps_lenient.log

if test -f ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps_lenient.g.vcf.gz
then
  echo "Lenient Hard Filter Complete: ./data/hard_filtered_variants/${sample_name}_hard_filtered_snps_lenient.g.vcf.gz $(date +%F)$" >> \
  ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_snps_lenient.log
fi

print "fin"
#################################

# Seperate samples

# Files:
# Snps
#./data/hard_filtered_variants/germline_hard_filtered_snps_lenient.g.vep.vcf
#./data/hard_filtered_variants/002683_Line-6_hard_filtered_snps_lenient.g.vep.vcf
#./data/hard_filtered_variants/002684_Line-7_hard_filtered_snps_lenient.g.vep.vcf
# Indels
#./data/hard_filtered_variants/germline_hard_filtered_indels_lenient.g.vep.vcf
#./data/hard_filtered_variants/002683_Line-6_hard_filtered_indels_lenient.g.vep.vcf
#./data/hard_filtered_variants/002684_Line-7_hard_filtered_indels_lenient.g.vep.vcf

# Seperate samples from gvcf file into respective vcf files
# Collective germline samples
find ./data/raw_snps_indels/ -name "[0-9]*-0*S*.g.vcf.gz" | \
xargs -i basename {} | \
sed 's/.g.vcf.gz//' | \
sort | uniq | \
xargs -i echo 'qsub ./scripts/seperate_samples.sh -v Var='{} |sh

# ./scripts/seperate_samples.sh
##############################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=00:30:00,mem=8gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# Select individual samples and filter out filtered calls and all non-variants
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/hard_filtered_variants/germline_hard_filtered_snps_lenient.g.vcf.gz \
-o ./data/hard_filtered_variants/${Var}_hard_filtered_snps_lenient.vcf.gz \
-sn ${Var} \
-env \
-ef

# env: Dont include non-variant sites
# ef: Don't include filtered sites
##############################

# Seperate samples from gvcf file into respective vcf files
# SNV calls on tumors (germline and somatic)
find /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/GATK_Variant_Analysis_Germline_DNA/data/raw_snps_indels/ -name "[7-9]*-[1-3]*S*raw_snps_indels.g.vcf" | head -n 25 | \
xargs -i basename {} | \
sed 's/_raw_snps_indels.g.vcf//' | \
sort | uniq | \
xargs -i echo 'qsub ./scripts/seperate_tumor_samples_SNV_calls.sh -v Var='{} |sh

# ./scripts/seperate_tumor_samples_SNV_calls.sh
##############################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=00:30:00,mem=8gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/GATK_Variant_Analysis_Germline_DNA/

# Select individual samples and filter out filtered calls and all non-variants
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/hard_filtered_variants/collective_samples_hard_filtered_snps_lenient.g.vcf \
-o ./data/hard_filtered_variants/${Var}_hard_filtered_snps_lenient.vcf.gz \
-sn ${Var} \
-env \
-ef

# env: Dont include non-variant sites
# ef: Don't include filtered sites
##############################

#/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/GATK_Variant_Analysis_Germline_DNA/data/hard_filtered_variants/collective_samples_hard_filtered_snps_lenient.g.vcf

# Seperate samples from gvcf file into respective vcf files
# lines 6 and 7
find ./data/hard_filtered_variants/ -name "00268[3-4]_Line-[6-7]_hard_filtered_snps_lenient.g.vcf.gz" | \
xargs -i basename {} | \
sed 's/_hard_filtered_snps_lenient.g.vcf.gz//' | \
sort | uniq | \
xargs -i echo 'qsub ./scripts/seperate_samples_parents.sh -v Var='{} |sh

# ./scripts/seperate_samples_parents.sh
##############################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=04:00:00,mem=8gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# Select individual samples and filter out filtered calls and all non-variants
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/hard_filtered_variants/${Var}_hard_filtered_snps_lenient.g.vcf.gz \
-o ./data/hard_filtered_variants/${Var}_hard_filtered_snps_lenient.vcf.gz \
-sn ${Var} \
-env \
-ef

# env: Dont include non-variant sites
# ef: Don't include filtered sites
##############################

# 1. Extract the indels from the call set

# line 6x7 F1
find ./data/raw_snps_indels -name "germline_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/extract_indels_from_genotyper.sh -v Var='{} |sh

# line 6
find ./data/raw_snps_indels -name "002683_Line-6_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/extract_indels_from_genotyper.sh -v Var='{} |sh

# line 7
find ./data/raw_snps_indels -name "002684_Line-7_raw_snps_indels_genotyped.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/extract_indels_from_genotyper.sh -v Var='{} |sh

# ./scripts/extract_indels_from_genotyper.sh
################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=03:59:00,mem=4gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

sample_name=$(basename ${Var} "_raw_snps_indels_genotyped.g.vcf.gz")

# 1. Extract the SNPs from the call set
java -Xmx4g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ${Var} \
-selectType INDEL \
-o ./data/raw_snps_indels/${sample_name}_raw_indels_extracted.g.vcf.gz \
2> ./analysis/raw_snps_indels/${sample_name}_raw_indels_extracted.log

if test -f ./data/raw_snps_indels/${sample_name}_raw_indels_extracted.g.vcf.gz
then
  echo "Raw SNP Extraction Complete: ./data/raw_snps_indels/${sample_name}_raw_indels_extracted.g.vcf.gz $(date +%F)" >> \
  ./analysis/raw_snps_indels/${sample_name}_raw_indels_extracted.log
fi

print "fin"
#################################

# 2. Determine parameters for filtering indels
# Extract QD Stats for R
python2.7 ./scripts/Collect_QD.py \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted.g.vcf \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted_QD.txt

# Extract FS Stats for R
python2.7 ./scripts/Collect_FS.py \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted.g.vcf \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted_FS.txt

# Extract MQ Stats for R
python2.7 ./scripts/Collect_MQ.py \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted.g.vcf \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted_MQ.txt

# Extract MQRankSum Stats for R
python2.7 ./scripts/Collect_MQRankSum.py \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted.g.vcf \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted_MQRankSum.txt

# Extract ReadPosRankSum Stats for R
python2.7 ./scripts/Collect_ReadPosRankSum.py \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted.g.vcf \
./data/hard_filtered_variants/collective_samples_raw_indels_extracted_ReadPosRankSum.txt


# Generate genomic statistic density plots in R
module load R/3.2.0
Rscript --vanilla ./scripts/hard_filt_plots_indels.R

# ./scripts/hard_filt_plots_indels.R
##########################################
# library("ggplot2")
# 
# # In R, Generate the QD Plot
# samples_QD = read.table("./data/hard_filtered_variants/collective_samples_raw_indels_extracted_QD.txt")
# qplot(V1, data = samples_QD, geom = "density") + xlab("QD") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_QD_density_indels.png")
# 
# # In R, Generate the FS Plot (x-axis is logged)
# samples_FS = read.table("./data/hard_filtered_variants/collective_samples_raw_indels_extracted_FS.txt")
# qplot(V1, data = samples_FS, geom = "density", log="x") + xlab("FS") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_FS_density_indels.png")
# # Cutoff right after right peak
# 
# # In R, Generate the MQ Plot
# samples_MQ = read.table("./data/hard_filtered_variants/collective_samples_raw_indels_extracted_MQ.txt")
# qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_indels.png")
# 
# # Zoom in on the MQ peak at 40
# qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(39.0, 41.0)
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_40_indels.png")
# 
# # Zoom in on the MQ peak at 60
# qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(59.0, 61.0)
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_60_indels.png")
# # Discard anything that is not 60
# 
# # In R, Generate the MQRankSum Plot
# samples_MQRankSum = read.table("./data/hard_filtered_variants/collective_samples_raw_indels_extracted_MQRankSum.txt")
# qplot(V1, data = samples_MQRankSum, geom = "density") + xlab("MQRankSum") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQRankSum_density_indels.png")
# # Remove anything less than -2
# 
# # In R, Generate the ReadPosRankSum Plot
# samples_ReadPosRankSum = read.table("./data/hard_filtered_variants/collective_samples_raw_indels_extracted_ReadPosRankSum.txt")
# qplot(V1, data = samples_ReadPosRankSum, geom = "density") + xlab("ReadPosRankSum") + ylab("Density")
# ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_ReadPosRankSum_density_indels.png")
##########################################


# 3. Apply the filter to the indel call set
# https://www.broadinstitute.org/gatk/guide/article?id=6925
# https://www.broadinstitute.org/gatk/guide/article?id=2806

# line 6x7 F1
find ./data/raw_snps_indels -name "germline_raw_indels_extracted.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/hard_filter_indels.sh -v Var='{} |sh

# line 6
find ./data/raw_snps_indels -name "002683_Line-6_raw_indels_extracted.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/hard_filter_indels.sh -v Var='{} |sh

# line 7
find ./data/raw_snps_indels -name "002684_Line-7_raw_indels_extracted.g.vcf.gz" |\
xargs -i echo 'qsub ./scripts/hard_filter_indels.sh -v Var='{} |sh

# ./scripts/hard_filter_indels.sh
################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=03:59:00,mem=8gb
#PBS -j oe
set -e
set -u
set -o pipefail

module load GATK/3.5.0

cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

sample_name=$(basename ${Var} "_raw_indels_extracted.g.vcf.gz")

# Perform stringent hard filtering
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/ref/galgal5.fa \
-V ${Var} \
--filterExpression "MQ < 60.0 || MQ > 60.0 || QD < 17.0 || FS > 0.903 || MQRankSum < -1.0" \
--filterName "INDEL_Hard_Filter" \
-o ./data/hard_filtered_variants/${sample_name}_hard_filtered_indels_stringent.g.vcf.gz \
2> ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_indels_stringent.log

if test -f ./data/hard_filtered_variants/${sample_name}_hard_filtered_indels_stringent.g.vcf.gz
then
  echo "Stringent Hard Filter Complete: ./data/hard_filtered_variants/${sample_name}_raw_indels_extracted.g.vcf.gz $(date +%F)$" >> \
  ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_indels_stringent.log
fi

# Perform lenient hard filtering
java -Xmx8g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/ref/galgal5.fa \
-V ${Var} \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filterName "INDEL_Hard_Filter" \
-o ./data/hard_filtered_variants/${sample_name}_hard_filtered_indels_lenient.g.vcf.gz \
2> ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_indels_lenient.log

if test -f ./data/hard_filtered_variants/${sample_name}_hard_filtered_indels_lenient.g.vcf.gz
then
  echo "Lenient Hard Filter Complete: ./data/hard_filtered_variants/${sample_name}_raw_indels_extracted.g.vcf.gz $(date +%F)$" >> \
  ./analysis/hard_filtered_variants/${sample_name}_hard_filtered_indels_lenient.log
fi

print "fin"
#################################

#

# Hard filtering over
###########################################################################################################


# Now it is time to create Truth files which will be used for machine learning to call Germline variants.
# First we will focus on Germline SNPs.
# We will first look at a 600K affymetrix array that called germline SNPs and genotypes on 6 6x7 F1 birds. We will take only the SNPs/genotypes
# that are identical between all 6 birds. These 6 birds were not used in our experiment but are the same F1 cross from 2 highly inbred breeds;
# Therefore the thought is that these SNPs/genotypes should also be present in our samples. Of course, this will
# not be true 100% of the time, therefore, we will check which of the supposed Germline SNPs/genotypes appear consistantly in each sample
# and if and only if they are consistant will they be passed on to construct our truth data.

python2.7 \
./scripts/Extract_6x7F1_Calls_600K_MA.py \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Microarray.txt \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Consistant.txt

# ./scripts/Extract_6x7F1_Calls_600K_MA.py
#######################################
import sys

infile=sys.argv[1]

outfile=open(sys.argv[2], 'w')

for line in open(infile):
  parts = line.split('\t')
  if parts[8] == parts[9] == parts[10] == parts[11] == parts[12] == parts[13]:
    outfile.write(parts[1] +'\t'+ parts[2] +'\t'+ parts[8] + '\n')
print 'Finished'
#######################################

# Liftover the Fixed MA Germline SNPs from Galgal4 to Galgal5:

# First convert the file to a new file that is compatable with galgal4 to galgal5 liftover:
python2.7 \
./scripts/prep_liftover_600K_Fixed_SNPs_4to5.py \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Consistant.txt \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Fixed_RD5LO.txt

# ./scripts/prep_liftover_600K_Fixed_SNPs_4to5.py
#######################################
import sys

infile=sys.argv[1]

outfile=open(sys.argv[2], 'w')

for line in open(infile):
  if not line[0] == 'C':
    parts = line.split('\t')
    try:
      start = int(parts[1]) - 1
      new_start = str(start)
      outfile.write(parts[0] +'\t'+ new_start +'\t'+ parts[1] +'\t'+ parts[2]) # no need to insert line break as parts[2] already contains it
    except ValueError:
      pass
print 'Finished'
#######################################

# LiftOver from GG4 to GG5
./scripts/liftover_gg4_to_gg5.sh \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Fixed_RD5LO.txt \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Fixed_Galgal5.bed \
./data/truth_files/600K_Germline_SNPs_6x7_F1_Fixed_Galgal4to5unmapped.bed

# ./scripts/liftover_gg4_to_gg5.sh
################################
set -e
set -u
set -o pipefail

module load BEDTools/2.24.0

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

sample_name=$(basename $1 "_rd5lo.txt")

# Sort the files, really is removing extra empty lines and cleaning up file
bedtools sort -i \
$1 > \
./data/truth_files/600K_germline_snps_int.bed

# LiftOver from GG4 to GG5
./scripts/liftOver \
./data/truth_files/600K_germline_snps_int.bed \
./scripts/galgal4togalgal5.liftover.chn \
$2 \
$3 \
2> ./analysis/liftover/liftover_GG4_to_GG5_${sample_name}.log

rm ./data/truth_files/600K_germline_snps_int.bed

echo "fin"
################################

# Extract only the locations in the 600K microarray that have alternate alleles than the galgal5 reference genome
find `pwd` -name "600K_Germline_SNPs_6x7_F1_Fixed_Galgal5.bed" | \
xargs -i echo 'qsub ./scripts/Extract_600K_MA_True_SNPs_6x7.sh -v Var='{} |sh

# ./scripts/Extract_600K_MA_True_SNPs_6x7.sh
################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=6:00:00,mem=1gb
#PBS -j oe
set -e
set -u
set -o pipefail

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

python \
./scripts/extract_600K_ma_true_snps.py \
${Var} \
./data/ref/galgal5.fa \
./data/truth_files/600K_Germline_TRUE_SNPs_6x7_Galgal5.bed \
./data/truth_files/600K_Germline_NOT_SNPs_6x7_Galgal5.bed

echo "Output files:"
echo "True SNPs File: ./data/truth_files/600K_Germline_TRUE_SNPs_6x7_Galgal5.bed"
echo "Not SNPs File: ./data/truth_files/600K_Germline_NOT_SNPs_6x7_Galgal5.bed"

qstat -f ${PBS_JOBID}
################################

# ./scripts/extract_600K_ma_true_snps.py
################################
import os
import sys
import subprocess 

infile1=(sys.argv[1])
infile2=(sys.argv[2])

os.system("module load SAMTools/1.2")

for line_MA in open(infile1):
  parts_MA = line_MA.split('\t')
  Chr_MA = str(parts_MA[0])
  Pos_MA = str(parts_MA[2])
  SNP_MA1 = str(parts_MA[3][0])
  SNP_MA2 = str(parts_MA[3][1])
  Ref_Out = subprocess.Popen(["samtools faidx" + " " + infile2 + " " + Chr_MA + ":" + Pos_MA + "-" + Pos_MA], stdout=subprocess.PIPE, shell=True)
  (out, err) = Ref_Out.communicate()
  Ref_SNP = out.split('\n')[1]
  if SNP_MA1 != Ref_SNP or SNP_MA2 != Ref_SNP:
    with open(sys.argv[3], 'a') as outfile_SNPs:
      outfile_SNPs.write(line_MA)
      outfile_SNPs.close()
  elif SNP_MA1 == Ref_SNP and SNP_MA2 == Ref_SNP:
    with open(sys.argv[4], 'a') as outfile_Not_SNPs:
      outfile_Not_SNPs.write(line_MA)
      outfile_Not_SNPs.close()
  else:
    print("You missed a possibility")
print "Fin"
################################


# Divide the 600K MicroArray File By Chromosome and position
# Galgal5 Info: http://www.ncbi.nlm.nih.gov/genome/?term=gallus
python2.7 \
./scripts/Divide_CHR_POS_MA.py \
./data/truth_files/600K_Germline_TRUE_SNPs_6x7_Galgal5.bed

# ./scripts/Divide_CHR_POS_MA.py
################################
import os,sys
from os import path

os.system("cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA")

path = './data/truth_files/divide_ma/'

infile1=sys.argv[1]

MA_CHR_Dict = {}
with open (infile1) as MA_FILE:
  for MA_LINE in MA_FILE:
    parts_MA = MA_LINE.split('\t')
    MA_CHR = str(parts_MA[0])
    MA_POS = int(parts_MA[2])
    key = MA_CHR
    val = MA_CHR
    MA_CHR_Dict[(key)] = val
for MA_SNP in open(infile1):
  parts_MA = MA_SNP.split('\t')
  MA_CHR = str(parts_MA[0])
  MA_POS = int(parts_MA[2])
  MA_ALT = str(parts_MA[3])
  if MA_CHR == '1' or MA_CHR == '2' or MA_CHR == '3' or MA_CHR == '4' or MA_CHR == '5' or MA_CHR == '6' or MA_CHR == '7' or MA_CHR == '8' or MA_CHR == '9' or MA_CHR == '10' or MA_CHR == '11' or MA_CHR == '12' or MA_CHR == '13' or MA_CHR == '14' or MA_CHR == '15' or MA_CHR == '16' or MA_CHR == '17' or MA_CHR == '18' or MA_CHR == '19' or MA_CHR == '20' or MA_CHR == '21' or MA_CHR == '22' or MA_CHR == '23' or MA_CHR == '24' or MA_CHR == '25' or MA_CHR == '26' or MA_CHR == '27' or MA_CHR == '28' or MA_CHR == 'W' or MA_CHR == 'Z':
    outfile1=open(path + '600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(MA_CHR) + '_POS_50000000.bed', 'a')
    outfile2=open(path + '600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(MA_CHR) + '_POS_100000000.bed', 'a')
    outfile3=open(path + '600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(MA_CHR) + '_POS_150000000.bed', 'a')
    outfile4=open(path + '600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(MA_CHR) + '_POS_200000000.bed', 'a')
    if MA_POS >= 0 and MA_POS < 50000000:
      outfile1.write(MA_SNP)
    if MA_POS >= 50000000 and MA_POS < 100000000:
      outfile2.write(MA_SNP)
    if MA_POS >= 100000000 and MA_POS < 150000000:
      outfile3.write(MA_SNP)
    if MA_POS >= 150000000 and MA_POS < 200000000:
      outfile4.write(MA_SNP)
    outfile1.close()
    outfile2.close()
    outfile3.close()
    outfile4.close()
print ('Fin')
################################


qsub ./scripts/Submit_Divide_CHR_POS_GATK.sh

# ./scripts/Submit_Divide_CHR_POS_GATK.sh
#################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=4gb
#PBS -j oe

module load Python/2.7.2

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

python2.7 \
./scripts/Divide_CHR_POS_GATK.py \
./data/hard_filtered_variants/collective_samples_raw_snps_extracted.g.vcf

qstat -f ${PBS_JOBID}
#################################

# ./scripts/Divide_CHR_POS_GATK.py
#################################
import os,sys
from os import path

os.system("cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA")

path = './data/truth_files/divide_ma/'

infile1=sys.argv[1]

for DNA_SNP in open(infile1):
  if not DNA_SNP[0] == '#':
    parts_DNA = DNA_SNP.split('\t')
    DNA_CHR = str(parts_DNA[0])
    DNA_POS = int(parts_DNA[1])
    DNA_ALT = str(parts_DNA[4])
    if DNA_CHR == '1' or DNA_CHR == '2' or DNA_CHR == '3' or DNA_CHR == '4' or DNA_CHR == '5' or DNA_CHR == '6' or DNA_CHR == '7' or DNA_CHR == '8' or DNA_CHR == '9' or DNA_CHR == '10' or DNA_CHR == '11' or DNA_CHR == '12' or DNA_CHR == '13' or DNA_CHR == '14' or DNA_CHR == '15' or DNA_CHR == '16' or DNA_CHR == '17' or DNA_CHR == '18' or DNA_CHR == '19' or DNA_CHR == '20' or DNA_CHR == '21' or DNA_CHR == '22' or DNA_CHR == '23' or DNA_CHR == '24' or DNA_CHR == '25' or DNA_CHR == '26' or DNA_CHR == '27' or DNA_CHR == '28' or DNA_CHR == 'W' or DNA_CHR == 'Z':
      outfile1=open(path + 'collective_samples_snps_chr_' + str(DNA_CHR) + '_pos_50000000.vcf', 'a')
      outfile2=open(path + 'collective_samples_snps_chr_' + str(DNA_CHR) + '_pos_100000000.vcf', 'a')
      outfile3=open(path + 'collective_samples_snps_chr_' + str(DNA_CHR) + '_pos_150000000.vcf', 'a')
      outfile4=open(path + 'collective_samples_snps_chr_' + str(DNA_CHR) + '_pos_200000000.vcf', 'a')
      if DNA_POS >= 0 and DNA_POS < 50000000:
        outfile1.write(DNA_SNP)
      if DNA_POS >= 50000000 and DNA_POS < 100000000:
        outfile2.write(DNA_SNP)
      if DNA_POS >= 100000000 and DNA_POS < 150000000:
        outfile3.write(DNA_SNP)
      if DNA_POS >= 150000000 and DNA_POS < 200000000:
        outfile4.write(DNA_SNP)
      outfile1.close()
      outfile2.close()
      outfile3.close()
      outfile4.close()
print ('Fin')

qstat -f ${PBS_JOBID}
#################################


# Make all the qsub files for the many versions of Compare_SNPs_Bed_VCF.py
python2.7 \
./scripts/make_compare_snps_scripts.py

# ./scripts/make_compare_snps_scripts.py
########################################
import os,sys
from os import path

os.system("cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA")

path = './scripts/compare_snps/'

for Chr in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', 'W', 'Z'):
  outfile=open(path + 'Compare_SNPs_scripts_Chr_' + str(Chr) + '_50000000.sh', 'w')
  outfile.write(r'#!/bin/bash -login' + '\n')
  outfile.write(r'#PBS -l nodes=1:ppn=1,walltime=15:00:00,mem=1gb' + '\n')
  outfile.write(r'#PBS -j oe' + '\n' + '\n' + '\n')
  outfile.write(r'module load Python/2.7.2' + '\n' + '\n')
  outfile.write(r'cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA' + '\n' + '\n')
  outfile.write('python2.7 \\' + '\n')
  outfile.write('./scripts/Compare_SNPs_Bed_VCF.py \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(Chr) + '_POS_50000000.bed \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_snps_chr_' + str(Chr) + '_pos_50000000.vcf \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_600K_SHARED_SNPs_CHR_' + str(Chr) + '_POS_50000000.vcf' + '\n' + '\n')
  outfile.write('qstat -f ${PBS_JOBID}' + '\n' + '\n')
  outfile.close()
  outfile=open(path + 'Compare_SNPs_scripts_Chr_' + str(Chr) + '_100000000.sh', 'w')
  outfile.write(r'#!/bin/bash -login' + '\n')
  outfile.write(r'#PBS -l nodes=1:ppn=1,walltime=15:00:00,mem=1gb' + '\n')
  outfile.write(r'#PBS -j oe' + '\n' + '\n' + '\n')
  outfile.write(r'module load Python/2.7.2' + '\n' + '\n')
  outfile.write(r'cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA' + '\n' + '\n')
  outfile.write('python2.7 \\' + '\n')
  outfile.write('./scripts/Compare_SNPs_Bed_VCF.py \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(Chr) + '_POS_100000000.bed \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_snps_chr_' + str(Chr) + '_pos_100000000.vcf \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_600K_SHARED_SNPs_CHR_' + str(Chr) + '_POS_100000000.vcf' + '\n' + '\n')
  outfile.write('qstat -f ${PBS_JOBID}' + '\n' + '\n')
  outfile.close()
  outfile=open(path + 'Compare_SNPs_scripts_Chr_' + str(Chr) + '_150000000.sh', 'w')
  outfile.write(r'#!/bin/bash -login' + '\n')
  outfile.write(r'#PBS -l nodes=1:ppn=1,walltime=15:00:00,mem=1gb' + '\n')
  outfile.write(r'#PBS -j oe' + '\n' + '\n' + '\n')
  outfile.write(r'module load Python/2.7.2' + '\n' + '\n')
  outfile.write(r'cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA' + '\n' + '\n')
  outfile.write('python2.7 \\' + '\n')
  outfile.write('./scripts/Compare_SNPs_Bed_VCF.py \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(Chr) + '_POS_150000000.bed \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_snps_chr_' + str(Chr) + '_pos_150000000.vcf \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_600K_SHARED_SNPs_CHR_' + str(Chr) + '_POS_150000000.vcf' + '\n' + '\n')
  outfile.write('qstat -f ${PBS_JOBID}' + '\n' + '\n')
  outfile.close()
  outfile=open(path + 'Compare_SNPs_scripts_Chr_' + str(Chr) + '_200000000.sh', 'w')
  outfile.write(r'#!/bin/bash -login' + '\n')
  outfile.write(r'#PBS -l nodes=1:ppn=1,walltime=15:00:00,mem=1gb' + '\n')
  outfile.write(r'#PBS -j oe' + '\n' + '\n' + '\n')
  outfile.write(r'module load Python/2.7.2' + '\n' + '\n')
  outfile.write(r'cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA' + '\n' + '\n')
  outfile.write('python2.7 \\' + '\n')
  outfile.write('./scripts/Compare_SNPs_Bed_VCF.py \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_' + str(Chr) + '_POS_200000000.bed \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_snps_chr_' + str(Chr) + '_pos_200000000.vcf \\' + '\n')
  outfile.write('./data/truth_files/divide_ma/collective_samples_600K_SHARED_SNPs_CHR_' + str(Chr) + '_POS_200000000.vcf' + '\n' + '\n')
  outfile.write('qstat -f ${PBS_JOBID}' + '\n' + '\n')
  outfile.close()
print "fin"
########################################


# Submit all the versions of Compare_SNPs_Bed_VCF.sh
find `pwd` -name "Compare_SNPs_scripts_Chr_1_*.sh" | \
xargs -i echo 'qsub '{} |sh

# ./scripts/compare_snps/Compare_SNPs_scripts_Chr_1_50000000.sh
########################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=15:00:00,mem=1gb
#PBS -j oe

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

module load Python/2.7.2

# Compare all of the SNPs separetely based on their chromosome and position
python2.7 \
./scripts/Compare_SNPs_Bed_VCF.py \
./data/truth_files/divide_ma/600K_Germline_SNPs_6x7_F1_Fixed_Galgal5_CHR_1_POS_50000000.bed \
./data/truth_files/divide_ma/collective_samples_snps_chr_1_pos_50000000.vcf \
./data/truth_files/divide_ma/collective_samples_600K_SHARED_SNPs_CHR_1_POS_50000000.vcf

qstat -f ${PBS_JOBID}
########################################

# ./scripts/Compare_SNPs_Bed_VCF.py
########################################
import sys

infile1=sys.argv[1]
infile2=sys.argv[2]

outfile=open(sys.argv[3], 'w')

for DNA_SNP in open(infile2):
        if not DNA_SNP[0] == '#':
                parts_DNA = DNA_SNP.split('\t')
                DNA_CHR_POS = (parts_DNA[0] + parts_DNA[1])
                DNA_ALT = parts_DNA[4]
                for MA_SNP in open(infile1):
                        parts_MA = MA_SNP.split('\t')
                        MA_CHR_POS = (parts_MA[0] + parts_MA[2])
                        MA_ALT = parts_MA[3]
                        MA_ALT = MA_ALT.replace("\n", "")
                        if MA_CHR_POS == DNA_CHR_POS and (DNA_ALT == MA_ALT[0] or DNA_ALT == MA_ALT[1]):
                                print (DNA_SNP)
                                outfile.write(DNA_SNP)
print('Fin')
########################################


grep "^#" ./data/hard_filtered_variants/collective_samples_raw_snps_extracted.g.vcf > \
./data/truth_files/collective_samples_600K_shared_snps.g.vcf

find -name "collective_samples_600K_SHARED_SNPs_CHR_*.vcf" | xargs cat >> ./data/truth_files/collective_samples_600K_shared_snps.g.vcf

module load vcftools/0.1.12a
module load tabix/0.2.6

bgzip ./data/truth_files/collective_samples_600K_shared_snps.g.vcf
vcf-sort ./data/truth_files/collective_samples_600K_shared_snps.g.vcf.gz > ./data/truth_files/collective_samples_600K_shared_snps.g.vcf
bgzip -c ./data/truth_files/collective_samples_600K_shared_snps.g.vcf > ./data/truth_files/collective_samples_600K_shared_snps.g.vcf.gz
tabix -p vcf ./data/truth_files/collective_samples_600K_shared_snps.g.vcf.gz

# Truth value dataset complete!
# SNPs found independently in DNASeq and Affymetrix microarrays found:
# ./data/truth_files/collective_samples_600K_shared_snps.g.vcf
# ./data/truth_files/collective_samples_600K_shared_snps.g.vcf.gz

###########################################################################################################


# We have generated a list of SNPs that were called independently with two seperate technologies 
# (although with different samples): Affymetrix SNP microarrays and DNAsequencing. We will treat this list
# as SNPs that we assume are true. We will then use these SNPs in the VSQR machine learning approach (Gaussian Mixture Model)
# to better generate a high confidence germline SNP dataset.

# Overview website on VSQR
# https://www.broadinstitute.org/gatk/guide/article?id=39

# Step-by-step documentation on performing VSQR
#https://www.broadinstitute.org/gatk/guide/article?id=2805

# FAQ on VSQR
# https://www.broadinstitute.org/gatk/guide/article?id=1259

#Steps
# Prepare recalibration parameters for SNPs
# a. Specify which call sets the program should use as resources to build the recalibration model
# b. Specify which annotations the program should use to evaluate the likelihood of Indels being real
# c. Specify the desired truth sensitivity threshold values that the program should use to generate tranches
# d. Determine additional model parameters
# Build the SNP recalibration model
# Apply the desired level of recalibration to the SNPs in the call set

# 1. Prepare recalibration parameters for SNPs

# a. Specify which call sets the program should use as resources to build the recalibration model

# True sites training resource: Shared SNPs between 600K MA
# /data/truth_files/collective_samples_600K_shared_snps_formated.g.vcf
# truth=true (representative of true sites)
# training=true (will be used to train the recalibration model)
# Likelihood of these variants should be around 97% or Q15 (96.84%) (Phred scale)

# Known sites resource, not used in training: dbSNP file
# ./data/vsqr/dbsnp.galgal5.vcf
# truth=false (Has not been validated to a high degree of confidence)
# training=false (Will not use for training)
# known=true (Stratifies output metrics like Ti/Tv ratio by whether variants are present in dbSNP or not)
# Likelihood=Q2 (36.90%)

# Transfer dbSNP file that was previously lifted over from Galgal4 to Galgal5 from TUM cluster to MSU cluster
# Files transferred from TUM cluster to MSU cluster on Mon Apr 25 16:41:02 CEST 2016
# rsync -avp --bwlimit=50000 /home/proj/MDW_genomics/xu/dbSNP/* \
# steepale@rsync.hpcc.msu.edu:/mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA/data/vsqr/

# b. Specify which annotations the program should use to evaluate the likelihood of SNPs being real

# https://www.broadinstitute.org/gatk/guide/article?id=1259
# The InbreedingCoeff is a population level statistic that requires at least 10 samples in order to be computed. 
# For projects with fewer samples, or that includes many closely related samples (such as a family) please omit this 
# annotation from the command line.

# Our samples are all from the same family, therefore, we will omit the InbreedingCoeff.

# c. Specify the desired truth sensitivity threshold values that the program should use to generate tranches
# We will calculate 4 tranches: 100.0, 99.9, 99.0, 90.0
# Based on graph(s) on this page: https://www.broadinstitute.org/gatk/guide/article?id=39

# Reformat resource files appropriately
# format of ./data/truth_files/collective_samples_600K_shared_snps.g.vcf
# grep -A1 "^#CHROM" ./data/truth_files/collective_samples_600K_shared_snps.g.vcf | column -t
# #CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  738-0_S33 741-0_S34 756-0_S35 766-0_S36 777-0_S37 
# 787-0_S38 788-0_S39 794-0_S40 798-0_S41 833-0_S42 834-0_S43 835-0_S44 841-0_S45 842-0_2_S28 
# 855-0_S46 863-0_S47 884-0_S48 901-0_2_S29 906-0_S49 911-0_S30 918-0_S50 927-0_S51
# 1 27599 . T C 730.06  SNP_Hard_Filter 
# AC=10;AF=0.357;AN=28;BaseQRankSum=0.358;ClippingRankSum=-7.360e-01;DP=57;ExcessHet=2.3480;FS=0.000;InbreedingCoeff=-0.1416;MLEAC=11;MLEAF=0.393;MQ=35.70;MQRankSum=0.736;QD=21.47;ReadPosRankSum=0.736;SOR=0.836  
# GT:AD:DP:GQ:PL  ./.:0,0:0 ./.:0,0:0 ./.:0,0:0 ./.:1,0:1 0/0:4,0:4:12:0,12,162 ./.:2,0:2 
# 0/1:1,2:3:23:77,0,23  0/1:3,3:6:99:103,0,109  1/1:0,4:4:12:165,12,0 0/1:1,2:3:34:72,0,34  0/1:1,2:3:34:76,0,34  
# 0/0:3,0:3:9:0,9,110 ./.:2,0:2 0/0:2,0:2:6:0,6,86  ./.:1,0:1 0/1:4,4:8:99:145,0,132  0/0:3,0:3:0:0,0,36  
# 0/1:3,2:5:70:70,0,103 1/1:0,2:2:6:86,6,0  ./.:2,0:2 0/0:1,0:1:3:0,3,38  0/0:2,0:2:6:0,6,76

# grep -v "^#" ./data/truth_files/collective_samples_600K_shared_snps.g.vcf | awk -F "\t" '{print NF; exit}'
# 31

# Need to sort again with picard
# Picard version: 1.131
java -jar ./scripts/picard/picard.jar SortVcf \
INPUT=./data/truth_files/collective_samples_600K_shared_snps.g.vcf \
OUTPUT=./data/truth_files/collective_samples_600K_shared_snps_sorted.g.vcf

# Extract the header then cut out the first 8 columns and feed it into a new file
(grep "^##" ./data/truth_files/collective_samples_600K_shared_snps_sorted.g.vcf; \
grep -v "^##" ./data/truth_files/collective_samples_600K_shared_snps_sorted.g.vcf | \
cut -f 1-8) > ./data/truth_files/collective_samples_600K_shared_snps_formated.g.vcf

# 2. Build the SNP recalibration model
# Make sure not to use the inbreeding coefficient as a parameter because all of our samples are from the same family
# https://www.broadinstitute.org/gatk/guide/article?id=1259

find `pwd` -name "collective_samples_raw_snps_extracted.g.vcf" | \
xargs -i echo 'qsub ./scripts/vsqr.sh -v Var='{} |sh
# AND
find `pwd` -name "S1-22_raw_snps_extracted.g.vcf" | \
xargs -i echo 'qsub ./scripts/vsqr.sh -v Var='{} |sh

# ./scripts/vsqr.sh
########################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=15:00:00,mem=60gb
#PBS -j oe

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

sample_name=$(basename ${Var} "_raw_snps_extracted.g.vcf")

module load GATK/3.5.0

# Build recalibration model
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ./data/ref/galgal5.fa \
-input ${Var} \
-resource:Affy_600K_MA_Matching_SNPs,known=false,training=true,truth=true,prior=15.0 ./data/truth_files/collective_samples_600K_shared_snps_formated.g.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./data/vsqr/dbsnp.galgal5.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.recal \
-tranchesFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.tranches \
-rscriptFile ./data/vsqr/${sample_name}_germline_recalibrate_snp_plots.R

# Apply recalibration at different levels of confidence
# 100.0% Confidence
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ./data/ref/galgal5.fa \
-input ${Var} \
-mode SNP \
--ts_filter_level 100.0 \
-recalFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.recal \
-tranchesFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.tranches \
-o ./data/vsqr/${sample_name}_germline_vsqr_100.0.g.vcf

# Apply recalibration
# 99.9% Confidence
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ./data/ref/galgal5.fa \
-input ${Var} \
-mode SNP \
--ts_filter_level 99.9 \
-recalFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.recal \
-tranchesFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.tranches \
-o ./data/vsqr/${sample_name}_germline_vsqr_99.9.g.vcf

# Apply recalibration
# 99.0% Confidence
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ./data/ref/galgal5.fa \
-input ${Var} \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.recal \
-tranchesFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.tranches \
-o ./data/vsqr/${sample_name}_germline_vsqr_99.0.g.vcf

# Apply recalibration
# 90.0% Confidence
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ./data/ref/galgal5.fa \
-input ${Var} \
-mode SNP \
--ts_filter_level 90.0 \
-recalFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.recal \
-tranchesFile ./data/vsqr/${sample_name}_germline_recalibrate_snp.tranches \
-o ./data/vsqr/${sample_name}_germline_vsqr_90.0.g.vcf

qstat -f ${PBS_JOBID}
########################################

# We have now generated four files of differring degrees of confidence of germline snp calls:
# The Germline samples from the tumor tissues:
# ./data/vsqr/collective_samples_germline_vsqr_100.0.g.vcf
# ./data/vsqr/collective_samples_germline_vsqr_99.9.g.vcf
# ./data/vsqr/collective_samples_germline_vsqr_99.0.g.vcf
# ./data/vsqr/collective_samples_germline_vsqr_90.0.g.vcf

# The Germline samples from the normal tissues:
# ./data/vsqr/S1-22_germline_vsqr_100.0.g.vcf
# ./data/vsqr/S1-22_germline_vsqr_99.9.g.vcf
# ./data/vsqr/S1-22_germline_vsqr_99.0.g.vcf
# ./data/vsqr/S1-22_germline_vsqr_90.0.g.vcf

# We have two files that represent true SNPs and false SNPs that originate from our Affymetrix 600K microarray comparison
# with the raw GATK calls on SNPs:
# File 1, ./data/truth_files/collective_samples_600K_shared_snps_formated.g.vcf (True positives): This file represents SNPs that were 
# consistant between all 6 F1 samples in our microarray calls, differed from the reference in at least 1 allele, and contained the same 
# alleles as the GATK raw SNP calls on atleast one germline sample.
# File 2, ./data/truth_files/600K_Germline_NOT_SNPs_6x7_Galgal5.bed

#####################################

# Genotype Refinement
#####################################

# We need to refine our genotypes.
# germline normal tissues are about 93% similar with the VSQR_100.g.vcf files; however, when sample A normal is compared
# to sample A tumor and sample A normal is compared to sample b tumor similarities in germline snps are 91.5% and 
# 90.5% respectively. 
# scripts for these analyses: /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_SNP_comparison_tumor_normal_tissues/compare_germline_snps_calls_tumor_normal.sh

# Website: https://www.broadinstitute.org/gatk/guide/article?id=4723
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CalculateGenotypePosteriors.php

# Combining variants from different files into one

# Normal tissue samples with the line 6 and line 7 germline calls
# https://www.broadinstitute.org/gatk/guide/article?id=53
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ./data/ref/galgal5.fa \
-V:line6,vcf ./data/vsqr/002683_Line-6_germline_vsqr_99.9.g.vcf \
-V:line7,vcf ./data/vsqr/002684_Line-7_germline_vsqr_99.9.g.vcf \
-V:normal,vcf ./data/vsqr/S1-22_germline_vsqr_99.9.g.vcf \
-o ./data/vsqr/normals_line6_line7_germline_vsqr_99.9.g.vcf

# Step 1: Derive posterior probabilities of genotypes (only family priors, no population priors because we lack these resources)
# Normal tissue samples with the line 6 and line 7 germline calls
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T CalculateGenotypePosteriors \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/normals_line6_line7_germline_vsqr_99.9.g.vcf \
--skipPopulationPriors \
-ped ./data/snp_genotype_refinement/normal_samples.ped \
-o ./data/vsqr/normals_family_priors_germline_vsqr_99.9.g.vcf

# Step 2: Filter low quality genotypes
# Normal tissue samples with the line 6 and line 7 germline calls
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/normals_family_priors_germline_vsqr_99.9.g.vcf \
-G_filter "GQ < 20.0" \
-G_filterName lowGQ \
-o ./data/vsqr/normals_filtered_gq_germline_vsqr_99.9.g.vcf

# Tumor tissue samples germline calls
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/collective_samples_germline_vsqr_99.9.g.vcf \
-G_filter "GQ < 20.0" \
-G_filterName lowGQ \
-o ./data/vsqr/tumors_filtered_gq_germline_vsqr_99.9.g.vcf

# Step 3: Annotate possible de novo mutations
java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/normals_filtered_gq_germline_vsqr_99.9.g.vcf \
-A PossibleDeNovo \
-ped ./data/snp_genotype_refinement/normal_samples.ped \
-o ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf

# Extract germline SNPs per sample and compare tumors with normal tissue counterparts

# Extract germline SNPs from normal tissue samples, parental line 6, and parental line 7 
# (note that lines 6 and 7 are not the actual parents and are likely male birds but highly inbred lines)

# Create array files with sample names

#Line 6 and Line 7 samples
grep -v "^##" ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | \
head -n 1 | \
cut -f 10,11 | \
tr '\t' '\n' > \
./data/arrays/parental_sample_names.txt

# Normal tissue samples
grep -v "^##" ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | \
head -n 1 | \
cut -f 12-33 | \
tr '\t' '\n' > \
./data/arrays/normal_tissue_sample_names.txt

# Tumor tissue samples
grep -v "^##" ./data/vsqr/tumors_filtered_gq_germline_vsqr_99.9.g.vcf | \
head -n 1 | \
cut -f 10-35 | \
tr '\t' '\n' > \
./data/arrays/tumor_tissue_sample_names.txt

./scripts/indv_samples/extract_snps_per_sample_normals.sh
###################################
#!/bin/bash -login
#PBS -l nodes=1:ppn=1,walltime=00:10:00,mem=5gb
#PBS -j oe
set -e
set -u
set -o pipefail

cd /mnt/scratch/steepale/birdman/MDV_project/DNAseq/GATK_Variant_Analysis_Germline_DNA

# Filter out true variants from samples by adequate genotype quality (GQ>20) and additional filters explained above (by FILTER tag)
declare -a norm
norm=( `cat "./data/arrays/normal_tissue_sample_names.txt" `)

for SAMPLE in "${norm[@]}"
do
  java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ./data/ref/galgal5.fa \
  -V ./data/germline_snps/normals_recal_famprioirs_gfiltered_denovo_snpeff_germline_99.9.g.vcf \
  -o './data/germline_snps/indv_samples/normals'${SAMPLE}'_normals_recal_famprioirs_gfiltered_denovo_snpeff_germline_99.9_int.g.vcf' \
  -sn ${SAMPLE} \
  --excludeNonVariants \
  --excludeFiltered

  grep -v "^##" \
  ./data/germline_snps/indv_samples/normals${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_snpeff_germline_99.9_int.g.vcf | \
  grep -v "lowGQ" > \
  ./data/germline_snps/indv_samples/normals${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_snpeff_germline_99.9.g.vcf

  rm ./data/germline_snps/indv_samples/normals${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_snpeff_germline_99.9_int.g.vcf
done

qstat -f ${PBS_JOBID}
####################################

declare -a rents
eval rents=(`cat "./data/arrays/parental_sample_names.txt" `)
for SAMPLE in "${rents[1]}"
do
java -Xmx5g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf \
-o ./data/germline_snps/indv_samples/parents/"${SAMPLE}"_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_int.g.vcf \
-sn "${SAMPLE}" \
--excludeNonVariants \
--excludeFiltered

grep -v "^##" \
./data/germline_snps/indv_samples/parents/"${SAMPLE}"_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_int.g.vcf | \
grep -v "lowGQ" > \
./data/germline_snps/indv_samples/parents/"${SAMPLE}"_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf

rm ./data/germline_snps/indv_samples/parents/"${SAMPLE}"_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_int.g.vcf
done

declare -a norm
eval norm=(` cat "./data/arrays/normal_tissue_sample_names.txt" `)
for SAMPLE in "${norm[21]}"
do
java -Xmx5g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf \
-o ./data/germline_snps/indv_samples/normals/${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_int.g.vcf \
-sn ${SAMPLE} \
--excludeNonVariants \
--excludeFiltered

grep -v "^##" \
./data/germline_snps/indv_samples/normals/${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_int.g.vcf | \
grep -v "lowGQ" > \
./data/germline_snps/indv_samples/normals/${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf

rm ./data/germline_snps/indv_samples/normals/${SAMPLE}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_int.g.vcf
done

declare -a tum
eval tum=(` cat "./data/arrays/tumor_tissue_sample_names.txt" `)
for SAMPLE in "${tum[@]}"
do
java -Xmx5g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./data/ref/galgal5.fa \
-V ./data/vsqr/tumors_filtered_gq_germline_vsqr_99.9.g.vcf  \
-o ./data/germline_snps/indv_samples/tumors/${SAMPLE}_tumors_gfiltered_germline_99.9_int.g.vcf \
-sn ${SAMPLE} \
--excludeNonVariants \
--excludeFiltered

grep -v "^##" \
./data/germline_snps/indv_samples/tumors/${SAMPLE}_tumors_gfiltered_germline_99.9_int.g.vcf | \
grep -v "lowGQ" > \
./data/germline_snps/indv_samples/tumors/${SAMPLE}_tumors_gfiltered_germline_99.9.g.vcf

rm ./data/germline_snps/indv_samples/tumors/${SAMPLE}_tumors_gfiltered_germline_99.9_int.g.vcf
done

# Compare the amount of line 6 and 7 snps in every normal and tumor matched pair

# Comparing Genotypes in both parents and tumor-normal samples
declare -a norm
norm=(`cat ./data/arrays/normal_tissue_sample_names.txt`)
declare -a rent
rent=(`cat ./data/arrays/parental_sample_names.txt`)
for p in {0..1}
do
echo "${rent[$p]} SNP amount:"
Psnp_tot=(`grep -v "^#" ./data/germline_snps/indv_samples/parents/${rent[$p]}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | wc -l`)
echo "${Psnp_tot}"
for n in {0..21}
do
echo "${norm[$n]} SNP amount:"
Nsnp_tot=(`grep -v "^#" ./data/germline_snps/indv_samples/normals/${norm[$n]}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | wc -l`)
echo "${Nsnp_tot}"
echo "Similar SNPs ${rent[$p]} and ${norm[$n]}:"
PT_snp_sim=(`comm -12 <(grep -v "^#" ./data/germline_snps/indv_samples/parents/${rent[$p]}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | cut -f 1,2,4,5,10 | awk -F ':' '{print $1}' | sort ) <(grep -v "^#" ./data/germline_snps/indv_samples/normals/${norm[$n]}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | cut -f 1,2,4,5,10 | awk -F ':' '{print $1}' | sort ) | wc -l`)
echo "${PT_snp_sim}"
Nsnp_ppt=$(( ${PT_snp_sim} * 100 / ${Nsnp_tot} ))
echo "${Nsnp_ppt}% of ${norm[$n]} SNPs are ${rent[$p]} SNPs"
Psnp_ppt=$(( ${PT_snp_sim} * 100 / ${Psnp_tot} ))
echo "${Psnp_ppt}% of ${rent[$p]} SNPs are ${norm[$n]} SNPs"
done
done

declare -a tum
tum=(`cat ./data/arrays/tumor_tissue_sample_names.txt`)
declare -a rent
rent=(`cat ./data/arrays/parental_sample_names.txt`)
for p in {0..1}
do
echo "${rent[$p]} SNP amount:"
Psnp_tot=(`grep -v "^#" ./data/germline_snps/indv_samples/parents/${rent[$p]}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | wc -l`)
echo "${Psnp_tot}"
for t in {0..25}
do
echo "${tum[$t]} SNP amount:"
Tsnp_tot=(`grep -v "^#" ./data/germline_snps/indv_samples/tumors/${tum[$t]}_tumors_gfiltered_germline_99.9.g.vcf | wc -l`)
echo "${Tsnp_tot}"
echo "Similar SNPs ${rent[$p]} and ${tum[$t]}:"
PT_snp_sim=(`comm -12 <(grep -v "^#" ./data/germline_snps/indv_samples/parents/${rent[$p]}_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf | cut -f 1,2,4,5,10 | awk -F ':' '{print $1}' | sort ) <(grep -v "^#" ./data/germline_snps/indv_samples/tumors/${tum[$t]}_tumors_gfiltered_germline_99.9.g.vcf | cut -f 1,2,4,5,10 | awk -F ':' '{print $1}' | sort ) | wc -l`)
echo "${PT_snp_sim}"
Tsnp_ppt=$(( ${PT_snp_sim} * 100 / ${Tsnp_tot} ))
echo "${Tsnp_ppt}% of ${tum[$t]} SNPs are ${rent[$p]} SNPs"
Psnp_ppt=$(( ${PT_snp_sim} * 100 / ${Psnp_tot} ))
echo "${Psnp_ppt}% of ${rent[$p]} SNPs are ${tum[$t]} SNPs"
done
done



