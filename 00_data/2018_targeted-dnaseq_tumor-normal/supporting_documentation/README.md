# Validation of Somatic SNVs and Indels (primarily non-coding) in *Gallus gallus* Marek's Disease (MD) Lymphomas via Targeted DNA Sequencing

We aim to validate somatic snv and indel candidates called from whole genome sequencing in MD lymphomas in line 6x7 chickens infected with 1000 pfu of Marek's Disease Virus (MDV) JM102W at hatch. There were 34 samples tested including a cohort of gonadal tumors and matched normals collected from 2014; additional tumors collected in 2014 and 2017 from gonad, heart, spleen, and liver. Targeted high coverage sequencing was performed on each region of interest by The RTSF Genomics Core at Michigan State University on one lane of an Illumina MiSeq (Project ID: STE7140). DNA was provided as PCR amplicons prepared at the Avian Disease and Oncology Lab. Each biological sample (1 well in 96 well plate) represents multiple amplicons from a single sample. These data represent raw reads (2x150bp paired end) from sequencing and the supporting documentation.

## Candidate Variant Pipeline

Variants were originally called from 22 MD gonadal tumors and matched normals in the illustrated pipeline.

<img src="./supporting_documentation/nonsynonymous_somatic_snv_indel_calling_pipeline.jpg" alt="drawing" width="1000px"/>

Statistical modeling revealed the most important variable to predict true variants to be the number of callers agreeing on the variant. Variants called by multiple callers, whether snvs or indels, were the most likely to be validated. 

## Supporting Documentation

-General report of assay from Genomics Core  
`supporting_documentation/20181128_SeqProduction_Cheng.xlsx`

-Supporting Email upon delivery
`supporting_documentation/Email_Carr_MiSeq_20181128_Amplicon_PE250.pdf`

-Genomics Core submission form  
`supporting_documentation/AmpliconProject_SubmissionIn96WellPlate.xlsx`

-Candidate variants  
`supporting_documentation/noncoding_variants_validation.txt`  

-Sample annotation  
`supporting_documentation/tissue_collection_and_sampling_2014_2017.txt`

-Amplicons and primers  
`supporting_documentation/amplicons_and_primers.xlsx`  

-Amplicon sequencing notes
`supporting_documentation/Amplicon_Sequencing_Notes.pdf`  

## Sequencing Files

All raw FASTQ sequences are in `fastq_files`

	$ find fastq_files -name "*.fastq.gz"
	fastq_files/017787-0_S27_L001_R1_001.fastq.gz
    fastq_files/017756-2_S10_L001_R1_001.fastq.gz
    fastq_files/017756-0_S1_L001_R2_001.fastq.gz
    fastq_files/017777-1_S22_L001_R2_001.fastq.gz

## Data Integrity

FastQC reports: `fastqc_raw`

MD5 file path: `fastq_files/fastq_raw_reads_amplicon-miseq.md5`

To perform MD5 check, from directory `fastq_files`

	$ md5sum -c fastq_raw_reads_amplicon-miseq.md5
	017756-0_S1_L001_R1_001.fastq.gz: OK
    017756-0_S1_L001_R2_001.fastq.gz: OK
    017756-1_S5_L001_R1_001.fastq.gz: OK
    017756-1_S5_L001_R2_001.fastq.gz: OK