# Calling Germline SNPs and Indels from Whole Genome Sequencing DNA-Seq 

Germline SNPs and Indels were called on parental lines 6 and 7, uninfected F1 6x7, and tumor and matched normal tissue from MDV infected F1s. We called germline SNPs and Indels to remove them from somatic SNV and Indel calls and to ensure that matched normal tissue was indeed from the same bird as its requisite tumor. Two general strategies were used; hard filtering (SNPs and Indels) and GATK's machine learning algorithm VSQR trained on 600K microarray data (SNPs only).

## Analyses in pipeline:

### Dependencies Used:
- Python (v.2.7.2)
- GATK (v3.5.0)
- Tabix (v0.2.6)
- BCFtools (v1.2)
- R (v3.2.0)
- GGplot2
- BEDTools (v2.24.0)
- Samtools (v1.2)
- VCFTools (v.0.1.12a)

### Formatting of Annotation and Reference:
- Construction of Truth dataset from f1 6x7 600K SNP arrays SNPs in [600K SNP Arrays of F1 6x7](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2012_affy-snp-600K_6-7-6x7):
    - Calls to line 6x7 were extracted that were all identical between the six samples of 6x7
    - Text file with calls converted to a BED format for preparation with liftOver
    - SNP coordinates lifted over from Galgal4 to Galgal5
    - SNPs were further selected for difference to Galgal5 reference genome
        - Galgal5 genome queried with Samtools faidx
    - SNPs were quired against hard filter SNPs and identical SNPs were considered truth set
    - Truth set sorted and processed with VCFtools and Tabix
    - Likelihood of these variants should be around 97% or Q15 (96.84%) (Phred scale)
- Known sites resource: [dbSNP file](https://github.com/steepale/IKZF1_paper_code/tree/master/01_reference_prep/dbSNP_gg4-gg5-lift)

### Quality Control:
- Variant counts were performed with BCFtools

### Major Steps in Pipeline Analysis:
- Raw SNPs and Indels were called from the GATK haplotype caller
- Joint genotyping was performed on all matched normal samples from infected F1s with GATK.
- Genotyping was performed on line 6 and line 7 birds, separately as single samples with GATK.
- Two filtering strategies were then applied...

- Hard filtering variants:
    - SNPs and Indels were extracted
    - The distribution of several statistics (QD, FS, MQ, MQRankSum, & ReadPosRankSum) were examined for hard filtering in R
    - Filters were applied to SNPs and Indels, separately with GATK.
        - SNPs (stringent): MQ < 60.0 || MQ > 60.0 || QD < 16.0 || FS > 0.778 || MQRankSum < 0.5
        - SNPs (lenient): QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0
        - Indels (stringent): MQ < 60.0 || MQ > 60.0 || QD < 17.0 || FS > 0.903 || MQRankSum < -1.0
        - Indels (lenient): QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0

- Filtering SNPs with the GATK VSQR (Gaussian mixture model)
    - Recalibration parameters were prepared for SNPs
        - Call sets were chosen as resources to build the recalibration model
            - Truth dataset from F1 6x7 600K SNP arrays SNPs
            - dbSNP
        - Truth sensitivity threshold values were calculated as 4 tranches: 100.0, 99.9, 99.0, and 90.0
    - The SNP recalibration model (VSQR--Gaussian mixture model) was build
- Finally, Genotype refinement was performed
    - We have family genotype information, so we will derive posterior probabilities of genotypes from family priors
    - The posterior probabilities of genotypes were built from family priors
    - Low quality genotypes (less than 20, Phred-scaled) were removed
    - Possible DeNovo mutations were annotated

## Data/Scripts:
A verbose workflow documentation script (should be adjusted in the future into a well designed workflow): `GATK_Variant_Analysis_Germline_DNA_Main_Documentation.sh`

## Genomic Datasets:
Training data:
- [2015_affy-snp_tumors-6-7-6x7](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2015_affy-snp_tumors-6-7-6x7)
- [2012_affy-snp-600K_6-7-6x7](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2012_affy-snp-600K_6-7-6x7)  

Whole Genome Sequencing:
- [2011_wgs_lines67-parental](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2011_wgs_lines67-parental)
- [2014_wgs_6x7-F1](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2014_wgs_6x7-F1)
- [2015_wgs_tumor-normal-gonad](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2015_wgs_tumor-normal-gonad)
 
