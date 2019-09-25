# Calling Germline SNPs and Indels from Whole Genome Sequencing DNA-Seq 

Germline SNPs and Indels were called on parental lines 6 and 7, uninfected F1 6x7, and tumor and matched normal tissue from MDV infected F1s. We called germline SNPs and Indels to remove them from somatic SNV and Indel calls and to ensure that matched normal tissue was indeed from the same bird as its requisite tumor. Two general strategies were used; hard filtering and GATK's machine learning algorithm VSQR trained on 600K microarray data. Germline variants were called from alignment files from [whole genome sequencing](https://github.com/steepale/IKZF1_paper_code/tree/master/02_fastq_to_bam/whole_genome_sequencing).

## Analyses in pipeline:

### Dependencies Used:
- GATK (v3.5.0)
- Tabix (v0.2.6)
- BCFtools (v1.2)
- R (v3.2.0)
- GGplot2

### Formatting of Annotation and Reference:
- Construction of Truth Dataset from informative SNPs in [600K SNP Arrays of F1 6x7](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2012_affy-snp-600K_6-7-6x7)

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

- Filtering with the GATK VSQR (Gaussian mixture model)


## Data/Scripts:


## Genomic Datasets:

- [600K SNP Arrays of F1 6x7](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2012_affy-snp-600K_6-7-6x7)

 
