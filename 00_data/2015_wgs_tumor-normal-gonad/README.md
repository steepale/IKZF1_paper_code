# Whole Genome Sequencing of Marek's Disease Lymphoma Gonadal Tumors and Matched Normal Tissues Infected with Marek's Disease Virus Strain JM102W

We collected lymphomas seeded in the gonads of males and females of F1 line 7 X line 6 (6x7) infected at hatch with 1000 pfu of MDV (strain JM102W). Birds were euthanized and tumors collected if birds were moribund or eight weeks of age. Whole genomic DNA was extracted via the QIAamp DNA Blood Mini Kit. Whole genome sequencing (WGS) was performed in 3 different batches; 26 gonadal tumor samples from 22 tumors from 22 birds on 12/30/2014, 3 matched normal samples on 05/28/2015, and 19 matched normal samples on 09/25/2015. All samples underwent WGS 125bp paired-end reads via the Illumina TruSeq Nano DNA Library Preparation Kit on Illumina HiSeq  machines at the Michigan State University (MSU) RTSF Genomics Core. Project designed and implemented by the Hans Cheng Lab (Hans.Cheng@ars.usda.gov) and Alec Steep (alec.steep@gmail.com).

## Sequencing Files

All raw FASTQ sequences are in `data`

		$ find data -name "*.fastq.gz"
		data/017756-3_CGGCTATG-CCTATCC_L002_R1_001.fastq.gz
		data/017842-2_2_TCCGCGAA-GTACTGA_L006_R2_001.fastq.gz
		data/017842-2_TCCGCGAA-CAGGACG_L008_R2_001.fastq.gz
		data/017835-1_TCCGCGAA-AGGCGAA_L007_R2_001.fastq.gz

## File Nomenclature

`017834-2_TCCGCGAA-CCTATCC_L001_R1_001.fastq.gz`  
`BIRD#-TUMOR#_BARCODE-BARCODE_LANE#_READ_001.fastq.gz`

`017834`: Bird ID  
`2`: Tumor ID relative to bird (e.g. tumor 2 collected from bird 017834)  
`TCCGCGAA-CCTATCC`: Sample Barcode (Demultiplexing across lanes)  
`L001`: Lane 1  
`R1`: Read 1

#### File nomeclature for replicate samples
`017834-2_2_TCCGCGAA-GGCTCTG_L001_R1_001.fastq.gz`  
`BIRD#-TUMOR#_Replicate_BARCODE-BARCODE_LANE#_READ_001.fastq.gz`

`017834`: Bird ID  
`2`: Tumor ID relative to bird (e.g. tumor 2 collected from bird 017834)  
`_2`: Identifier of biological replicate from distant side of tumor  
`TCCGCGAA-GGCTCTG`: Sample Barcode (Demultiplexing across lanes)  
`L001`: Lane 1  
`R1`: Read 1

## Data Integrity

MD5 file path: `data/20150211_md_gonad-tumors_6x7F1.md5`

To perform an MD5 check, in directory `data`

		$ md5sum -c 20150211_md_gonad-tumors_6x7F1.md5
		./017738-1_CGGCTATG-TATAGCC_L001_R1_001.fastq.gz: OK
		./017738-1_CGGCTATG-TATAGCC_L001_R2_001.fastq.gz: OK
		./017738-1_CGGCTATG-TATAGCC_L002_R1_001.fastq.gz: OK
		./017738-1_CGGCTATG-TATAGCC_L002_R2_001.fastq.gz: OK

## Supporting Annotation

Supporting annotation of birds, tumors, samples, and extracted DNA found in `docs`

Email explanations of sequencing from the MSU RTSF Genomics Core

	Email_from_Kevin_Carr_Notification_DNASeq_Complete_2015_01_15.pdf
	Email_from_Kevin_Carr_Notification_RNASeq_DNASeq_Complete_2015_06_03.pdf
	Email_from_Kevin_Carr_Notification_DNASeq_Complete_2015_10_13.pdf

Sequencing summaries by lane and sample

	2015_02_13_RTSF_Output_DNA_S1-26.xlsx
	2015_05_28_RTSF_Output_DNA_S28-30_RNA_S1.3.5.7.9-11.13-15.17.19.22.23.26.xlsx
	2015_09_25_RTSF_Output_DNA_S33-51.xlsx

RTSF submission forms

	Steep_RTSFSubmissionForm_Tubes_DNASeq_10.20.14.xlsx
	RTSF_Submission_Form_RNASeq-Tumor-Controls_DNASeq-Controls_03-12-15.xlsx
	Steep_RTSFSubmissionForm_Tubes_DNASeq_09_08_15.xlsx

Decision making and testing of genomic integrity of extracted DNA

	DNASeq_Purity_Yeild_Size_10.16.14.xlsx

Explanation of samples, tissue, genomics material, and analyses

	Tumor_Sampling_Info.xlsx

