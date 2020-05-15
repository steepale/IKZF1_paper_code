  #===============================================================================
#
#         FILE: /Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/scripts/20200103_gatk-600K-affy-calls-MSUHPCC_steep.sh
#
#        USAGE: Line by line, documentation script
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for Germline variant
#                calls (SNPs) on preprocessed BAM files (600K loci from microarry)
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2020.01.03
#     REVISION:  
#===============================================================================

# Script is used to perform targeted SNP calling on 600K microarry sites, included sites without ALT alleles.
# Stats per site used to aid hard filtering and machine learning filters.

# Change to working directory
cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# Pull the docker image (must name the image as steepale_gatk.sif)
singularity build steepale_gatk.sif docker://steepale/gatk:3.5

# Perform this anlysis in a docker image
singularity shell steepale_gatk.sif


# Call Germline raw variants (SNPs and Indels) on DNA Sequencing samples
for bam in `find /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/final_bam -name "*-0_S*_Bwa_RG_dedupped_realigned.bam"`
do
echo ${bam}
sbatch --mem=20000 --export=bam=${bam},gb=20 ./scripts/20191214_GATK-600K-affy-calls-haplotypecaller_steep.sb 
done

# ./scripts/20191214_GATK-600K-affy-calls-haplotypecaller_steep.sb
################################
#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=general
#SBATCH --partition=general-long

########## Command Lines to Run ##########

# Change to working directory
cd /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/germline_snps_indels

# VARIABLES
###########################
# Cellect sample basename
sample_name=$(basename ${bam} "_Bwa_RG_dedupped_realigned.bam")
# output raw snp vcf
raw_snp_vcf="./data/raw_snps_indels/${sample_name}-raw-snps-indels-600K-regions_steep.g.vcf"
raw_snp_vcf_gz="./data/raw_snps_indels/${sample_name}-raw-snps-indels-600K-regions_steep.g.vcf.gz"
# directories
bam_dir="/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/final_bam"
ref="/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/galgal5/galgal5.fa"
snp_list="./data/600K-regions_steep.list"

# COMMANDS
###########################

# Notification of sample
echo ${sample_name}

# Call SNPs at 600K loci
singularity exec steepale_gatk.sif \
java -Xmx${gb}g -jar /opt/gatk3.5/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${ref} \
-I ${bam} \
-ERC BP_RESOLUTION \
--interval_padding 100 \
-L ${snp_list} \
-o ${raw_snp_vcf} \
2> ./analysis/raw_snps_indels/${sample_name}_GATK-600K-affy-calls-haplotypecaller.log

# Compress the vcf file
singularity exec steepale_gatk.sif \
bgzip -f ${raw_snp_vcf}
# Index the vcf file
singularity exec steepale_gatk.sif \
tabix -p vcf ${raw_snp_vcf_gz}

##################################






# Genotype these loci collectively across samples
java -jar /opt/gatk3.5/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R ${ref} \
--never_trim_vcf_format_field \
--includeNonVariantSites \
-V ./data/raw_snps_indels/833-0_S42-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/901-0_S29-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/884-0_S48-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/741-0_S34-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/863-0_S47-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/842-0_S28-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/834-0_S43-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/787-0_S38-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/927-0_S51-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/777-0_S37-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/756-0_S35-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/918-0_S50-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/835-0_S44-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/911-0_S30-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/841-0_S45-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/738-0_S33-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/788-0_S39-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/855-0_S46-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/794-0_S40-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/766-0_S36-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/906-0_S49-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-V ./data/raw_snps_indels/798-0_S41-raw-snps-indels-600K-regions_steep.g.vcf.gz \
-o ./data/raw_snps_indels/collective-genotypes-600K-regions_steep.g.vcf.gz

