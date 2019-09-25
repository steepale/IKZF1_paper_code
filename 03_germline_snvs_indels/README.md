# Processing and Mapping of RNA-Seq Data for Differential Gene Expression Analysis in Chicken and Marek's Disease Virus (MDV) Transcriptome 

RNA-Seq reads were processed and mapped to an appended Galgal5 genome with the MDV genome as an additional chromosome. Annotation was also adjusted to include both genomes. The rational was to capture the reads from integrated copies of the MDV genome. Although we could not consistently determine the integration sites of MDV genomes--due to long repeat regions flanking integration and short reads from our sequencing approach--we were able to examine RNA quantities from MDV in tumor tissue. However, it is important to acknowledge that we cannot differentiate reads mapping to the MDV genomes as resulting from either cytolytic or genome-integrated-latent MDV.

## Analyses in pipeline:

### Dependencies Used:
- Python (v2.7.2)
- Java (v1.8.0_31)
- FastQC (v0.11.3)
- RSeQC (v2.4)
- Trimmomatic (v0.33)
- TopHat (v2.1.0)
- BowTie2 (v2.2.6)
- STAR (v2.5.1b)
- BBTools (v38.22)
- cufflinks (v2.2.1)
- Samtools (v1.2)
- PySAM (v0.6)
- HTSeq (v0.6.1)

### Formatting of Annotation and Reference:
- The Marek's disease virus (MDV) reference genome, gff3 annotation, and transcriptome were obtained from [NCBI](https://www-ncbi-nlm-nih-gov.proxy1.cl.msu.edu/nuccore/NC_002229.3)
- The MDV annotation file is only available as a gff3 file. We used a custom script (`MDV_gff3_to_gtf.py`) to convert the gff3 format to gtf format.
- The Galgal5 reference genome, gtf annotation, and transcriptome were also obtained form Ensembl (release 92) (we did not use [the previously prepared reference genome](https://github.com/steepale/IKZF1_paper_code/tree/master/01_reference_prep))
- The MDV genome was added to the Galgal5 genome as an extra chromosome
- The MDV and Galgal5 annotations were also combined
- Reference genomes and transcriptomes were indexed via STAR.

### Quality Control:
- Chicken (Galgal4) specific rRNA was obtained from GenBank (accession KT445934)
- Strandedness of RNA sequencing was assessed via RSeQC ([RNASeq is stranded](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2015_rnaseq_tumor-normal-gonad)).

### Major Steps in Pipeline Analysis:
In all steps of this analysis, only paired reads were considered; therefore, if a read lost its pair during any processing step, then it was dismissed.
- Raw reads were assessed for quality via FastQC.
- Reads were trimmed of adaptors and low quality reads via Trimmomatic.
- Trimmed reads that mapped to chicken rRNA were removed via BBDuk
- Processed reads were assessed for quality via FastQC.
- Processed reads were mapped via STAR.
- The number of mapped reads were counted via Samtools.
- Alignment sam files were compressed to bam files with Samtools.
- Alignment (bam) files were sorted seperately by read name and position with Samtools.
- Reads were counted per gene/exon via HTSeq (PySAM is dependency).

## Data/Scripts:
- Primary documentation script to perform improved mapping with combo genome and improved aligner: `chicken_MDV_genome_main_documentation.txt`
- Original attempt to map and count reads. This is not the main script, rather used for processing reads: `rna_gene_expression_analysis_main_documentation.txt`  
- Python script to convert gff3 MDV file to gtf file format: `MDV_gff3_to_gtf.py`

## Genomic Dataset:
- [2015_rnaseq_tumor-normal-gonad](https://github.com/steepale/IKZF1_paper_code/tree/master/00_genomic_datasets/2015_rnaseq_tumor-normal-gonad)

 
