# Calling Germline SNPs and Indels from Whole Genome Sequencing DNA-Seq 

Germline SNPs and Indels were called on parental lines 6 and 7, uninfected F1 6x7, and tumor and matched normal tissue from MDV infected F1s. We called germline SNPs and Indels to remove them from somatic SNV and Indel calls and to ensure that matched normal tissue was indeed from the same bird as its requisite tumor. Two general strategies were used; hard filtering and GATK's machine learning algorithm VSQR trained on 600K microarray data. Germline variants were called from alignment files from [whole genome sequencing](https://github.com/steepale/IKZF1_paper_code/tree/master/02_fastq_to_bam/whole_genome_sequencing).

## Analyses in pipeline:

### Dependencies Used:
- GATK v(3.5.0)


### Formatting of Annotation and Reference:


### Quality Control:


### Major Steps in Pipeline Analysis:
- Raw SNPs and Indels were called from the GATK haplotype caller
- Joint genotyping was performed on x samples with 


## Data/Scripts:


## Genomic Datasets:


 
