# Processing and Mapping of RNA-Seq Data for Somatic and Germline Variant Calling

Reads were mapped according to the GATK best practices and STAR manual in preparation to call somatic and germline variants from RNA-Seq aligned data.

## Analyses in pipeline:

### Dependencies Used:
- FastQC (v0.11.4)
- Trimmomatic (v0.35)
- sickle (v1.33)
- STAR (v2.5.1b)
- SortMeRNA v(2.1)
- Seqtk (v1.2-r95-dirty)
- Picard-Tools (v1.141)
- GATK (v3.6)

### Formatting of Annotation and Reference:
- The major contigs of the Galgal5 genome were used as a reference as well as the supporting gtf file [explained here](https://github.com/steepale/IKZF1_paper_code/tree/master/01_reference_prep)
- The reference was indexed with STAR
- The Galgal5 dbSNP file was obtained from TODO


### Quality Control:


### Major Steps in Pipeline Analysis:
- Paired-end reads were interleaved with Seqtk
- Ribosomal rRNA was filtered out of interleaved reads with SortMeRNA
- Interleaved reads were separated to respective pairs with Seqtk
- Reads were trimmed of adaptors and low quality bases with Trimmomatic
- Reads were mapped to the reference genome with STAR (2-pass mapping)
- Aligned sam files were converted to bam format with Picard-Tools
- Read groups were added with Picard-Tools
- Bam files were sorted by coordinate with Picard-Tools
- Duplicate reads were removed by Picard-Tools
- Mapped reads with N's in their CIGAR strings were split and mapping qualities were reassigned with GATK
- Mapped reads were realigned around Indels with GATK
- Base quality scores were recalibrated with GATK
- Alignment files were merged across lanes with Picard Tools and underwent another round of processing:
    - Mapped reads with N's in their CIGAR strings were split and mapping qualities were reassigned
    - Mapped reads were realigned around Indels
    - Base quality scores were recalibrated
    - Final bam files were indexed with Picard-Tools

## Data Files:
- Sample/sequencing-run config file: `config.txt`  
- Perl script to run pipeline: `rnaseq.pl`  
- Script to submit jobs to cluster environment: `submit.pl`

## Genomic Dataset:
- [2015_rnaseq_tumor-normal-gonad](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2015_rnaseq_tumor-normal-gonad)

### Notes:
Two samples from uninfected birds (017824 & 017748) did not undergo this analysis.

Analysis was inspired by the Genome Analysis Toolkit (GATK) best practices pipeline as well as the STAR manual.

Analysis designed and implemented by [Hongen Xu](https://github.com/hongenxu).