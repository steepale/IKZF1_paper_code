  #===============================================================================
#
#         FILE: /Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/scripts/20200103_gatk-600K-affy-calls_steep.sh
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

# Perform this anlysis in a docker image
docker run -it \
-v /Users/Alec/Documents/Bioinformatics/MDV_Project/galgal5:/Users/Alec/Documents/Bioinformatics/MDV_Project/galgal5 \
-v /Volumes/Frishman_4TB/bams_wgs_2014/normal:/Volumes/Frishman_4TB/bams_wgs_2014/normal \
-v /Users/Alec/Documents/Bioinformatics/MDV_Project/galgal5:/Users/Alec/Documents/Bioinformatics/MDV_Project/galgal5 \
-v /Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data:/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data \
-v /Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/raw_snps_indels:/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/raw_snps_indels \
--rm \
--name gatk \
steepale/gatk:3.5

# Change to working directory
cd /Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels

# Need to update docker file to have tabix
sudo apt-get install tabix

# Set a few variables
bam_dir="/Volumes/Frishman_4TB/bams_wgs_2014/normal/"
ref="/Users/Alec/Documents/Bioinformatics/MDV_Project/galgal5/galgal5.fa"
snp_list="./data/600K-regions-sub_steep.list"

# We call SNPs on the 600K loci to feed into our machine learning model (combo of hard filters and machine learning)
for bam in `ls -1 ${bam_dir}*-0_S*_Bwa_RG_dedupped_realigned.bam`
do
# VARIABLES
###########################
# Cellect sample basename
sample=`basename ${bam} | sed 's/_Bwa_RG_dedupped_realigned.bam//'`
# output raw snp vcf
raw_snp_vcf="./data/raw_snps_indels/${sample}-raw-snps-indels-600K-regions_steep.g.vcf"
raw_snp_vcf_gz="./data/raw_snps_indels/${sample}-raw-snps-indels-600K-regions_steep.g.vcf.gz"

# COMMANDS
###########################

echo ${sample}
# Call SNPs at 600K loci
java -jar /opt/gatk3.5/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${ref} \
-I ${bam} \
-ERC BP_RESOLUTION \
-L ${snp_list} \
-o ${raw_snp_vcf}

# Compress the vcf file
bgzip -f ${raw_snp_vcf}
tabix -p vcf ${raw_snp_vcf_gz}
#############################
done

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

# Exit the docker image
exit

# Remove the docker image




# TODO: consider interval padding
#--interval_padding 100 \
