# Validation of Nonsynonymous Somatic Mutations in *Gallus gallus* Marek's Disease (MD) Lymphomas via Agrixplex Targeted DNA Sequencing

We aim to validate somatic nonsynonymous snv and indel candidates (131 high priority candidates and 6 -- 50 bp regions of interest) called from whole genome sequencing in MD lymphomas in line 6x7 chickens infected with 1000 pfu of Marek's Disease Virus (MDV) JM102W at hatch. There were 188 samples tested including a cohort of gonadal tumors and matched normals collected from 2014; additional tumors collected in 2014 and 2017 from gonad, heart, spleen, thymus, pancreas, liver, proventiculus, kidney, and bursa; line 6 and line 7 (male and female) germline samples; and MD cell lines MSB1, RP2, and RP19. Targeted high coverage sequencing was performed on each variant and region of interest by Agriplex Genomics using the PlexSeq^(TM)^ method and analyzed with PlexCall^(TM)^ software. Whole genomic DNA was provided to Agriplex in 2 -- 96 well plates. These data represent raw reads from Agrixplex Genomics sequencing and the supporting documentation.

## Candidate Variant Pipeline

Variants were originally called from 22 MD gonadal tumors and matched normals in the illustrated pipeline.

<img src="./supporting_documentation/nonsynonymous_somatic_snv_indel_calling_pipeline.jpg" alt="drawing" width="1000px"/>

Statistical modeling revealed the most important variable to predict true variants to be the number of callers agreeing on the variant. Variants called by 4 or more callers, whether snvs or indels, were the most likely to be validated. 

## Supporting Documentation

-General report of assay from Agriplex  
`supporting_documentation/AM-MSU_Cheng_Chicken_PlexSeq_Report_090517.docx`

-Agriplex submission form  
`supporting_documentation/SNP_Submission_form_top150_w50bp.xlsx`

-Candidate variants  
`supporting_documentation/somatic_snvs_and_indels_final_priorityv2.xlsx`  
`supporting_documentation/somatic_snvs_and_indels_final_priorityv2.txt`

-Sample annotation  
`supporting_documentation/MD_tumor_6x7_JM102W_annotation_database_2014_2017.xlsb`

-Sample and plate organization  
`supporting_documentation/Sample_Name_Submission_Platemap_Agriplex_v3.xlsb`

-Variant validation results from Agriplex  
`supporting_documentation/MSU_Cheng_Chicken_PlexSeq_Report_090517_SNPs_&_Indels.xlsx`  
`supporting_documentation/MSU_Cheng_Chicken_PlexSeq_Report_090517_Remaining_Indels.xlsx`  
`supporting_documentation/MSU_Cheng_Chicken_PlexSeq_Report_090517_SNPs_Indels.txt`

-Primer sequences flanking variants and regions of interest  
`supporting_documentation/100717_MSU_Chicken_PlexSeq_Primers.xlsx`  
`supporting_documentation/100717_MSU_Chicken_PlexSeq_Primers.txt`

## Sequencing Files

All raw FASTQ sequences are in `fastq_files`

	$ find fastq_files -name "*.fastq.gz"
	fastq_files/MSU_834-2-PL01-F03_S22_L001_R1_001.fastq.gz
	fastq_files/MSU_834-0-PL01-E03_S21_L001_R1_001.fastq.gz
	fastq_files/MSU_863-0-PL01-A05_S33_L001_R1_001.fastq.gz
	fastq_files/MSU_822-1-PL02-A07_S143_L001_R1_001.fastq.gz

## File Nomenclature

`MSU_101-T-PL01-C11_S83_L001_R1_001.fastq.gz`  
`MSU_BirdID-Tissue-Plate#-Well#_Sample#_Lane#_Read#_AdditionalID.fastq.gz`

`MSU`: Michigan State University "**Go Green!**"  
`101`: Bird 101  
`T`: Thymus tumor  
`PL01`: Plate 1  
`C11`: Well C11  
`S83`: Sample 83  
`L001`: Lane 1  
`R1`: Read 1  

## Data Integrity

FastQC reports: `fastqc_raw`

MD5 file path: `fastq_files/fastq_raw_reads_validation.md5`

To perform MD5 check, from directory `fastq_files`

	$ md5sum -c fastq_raw_reads_agrixplex_validation.md5
	MSU_101-T-PL01-C11_S83_L001_R1_001.fastq.gz: OK
	MSU_105-H-PL02-C02_S105_L001_R1_001.fastq.gz: OK
	MSU_108-P-PL02-D02_S106_L001_R1_001.fastq.gz: OK
	MSU_109-G-PL01-D09_S68_L001_R1_001.fastq.gz: OK


