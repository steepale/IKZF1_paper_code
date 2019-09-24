# Processing and Mapping of Whole Genome Sequencing Datasets

## Analyses in pipeline:
Reads were inspected for quantity and quality before and after trimming with FastQC (v0.11.4). Reads were trimmed of adaptors and low quality bases with Trimmomatic (v0.35) and sickle (v1.33), respectively. Trimmed reads were aligned to the Gallus gallus 5 reference genome (major contigs) with BWA-MEM (v0.7.12-r1044). Read group annotation was added via Picard tools (v1.141). Reads within each sample and sequencing lane were filtered of duplicate reads and were realigned around Indels via Picard tools and GATK (v3.5). Indel realignment was performed again after samples were merged across lanes.

## Data files:
- Sample/sequencing-run config file: `config.txt`  
- Perl script to run pipeline: `fastq2bam.pl`  
- Python(ic) script for improved pipeline (not used in IKZF1 paper): `fastq2bam_main_documentation.txt`
- Script to submit jobs to cluster environment: `submit.pl`

## Genomic datasets:
- [2011_wgs_lines67-parental](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2011_wgs_lines67-parental)
- [2014_wgs_6x7-F1](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2014_wgs_6x7-F1)
- [2015_wgs_tumor-normal-gonad](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2015_wgs_tumor-normal-gonad)


Analysis was inspired by the Genome Analysis Toolkit (GATK) best practices pipeline.

Perl script written for Technical University of Munich cluster by [Hongen Xu](https://github.com/hongenxu).
