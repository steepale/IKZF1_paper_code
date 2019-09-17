# RNA-Seq (mRNA) of Marek's Disease (MD) Lymphoma Gonadal Tumors and Isolated Uninfected CD4+ T-cells (similar to matched normal)

We collected MD lymphomas seeded in the gonads of males and females of F1 birds line 7 x line 6 (6x7) infected at hatch with 1000 pfu of MDV (strain JM102W). Birds were euthanized and tumors collected when birds were moribund or if they reached eight weeks of age (birds challenged at hatch). RNA was extracted via the miTNeasy mini kit (Qiagen). Sequencing was performed in 1 batch; 26 gonadal tumor samples from 22 tumors from 22 birds and 8 samples of isolated CD4+ spleenic T-cells from 8 uninfected birds on 05/28/2015. RNA-Seq libraries were prepared using the NuGen Ovation Single Cell RNA-Seq System (stranded). All samples underwent sequencing to produce 125bp paired-end reads via Illumina HiSeq  machines at the Michigan State University (MSU) RTSF Genomics Core. Project designed and implemented by the Hans Cheng Lab (Hans.Cheng@ars.usda.gov) and Alec Steep (alec.steep@gmail.com).

## Sequencing Files

All raw FASTQ sequences are in `data`

		$ find data -name "*.fastq.gz"
		data/017884-2_AACCAG_L003_R2_001.fastq.gz
		data/017777-3_AAGAGG_L001_R2_001.fastq.gz
		data/017855-1_1_CGTAGA_L003_R1_001.fastq.gz
		data/017834-2_GTCGTA_L005_R2_001.fastq.gz

## File Nomenclature

`017834-2_GTCGTA_L005_R2_001.fastq.gz`  
`BIRD#-TUMOR#_BARCODE_LANE#_READ_001.fastq.gz`

`017834`: Bird ID  
`2`: Tumor ID relative to bird (e.g. tumor 2 collected from bird 017834)  
`GTCGTA`: Sample Barcode (Demultiplexing across lanes)  
`L005`: Lane 5  
`R2`: Read 2

#### File nomeclature for replicate samples
`017901-2_2_TGGTGA_L001_R1_001.fastq.gz`  
`BIRD#-TUMOR#_Replicate_BARCODE-BARCODE_LANE#_READ_001.fastq.gz`

`017901`: Bird ID  
`2`: Tumor ID relative to bird (e.g. tumor 2 collected from bird 017901)  
`_2`: Identifier of biological replicate from distant side of tumor  
`TGGTGA`: Sample Barcode (Demultiplexing across lanes)  
`L001`: Lane 1  
`R1`: Read 1

## Data Integrity

MD5 file path: `data/20150528_md_gonad-tumors_6x7F1_RNASeq.md5`

To perform an MD5 check, in directory `data`

		$ md5sum -c 20150528_md_gonad-tumors_6x7F1_RNASeq.md5
		./017884-2_AACCAG_L003_R2_001.fastq.gz : OK
		./017777-3_AAGAGG_L001_R2_001.fastq.gz : OK
		./017855-1_1_CGTAGA_L003_R1_001.fastq.gz : OK
		./017834-2_GTCGTA_L005_R2_001.fastq.gz : OK

## Supporting Annotation

Supporting annotation of birds, tumors, samples, and extracted DNA found in `docs`

Email explanations of sequencing from the MSU RTSF Genomics Core

	Email_from_Kevin_Carr_Notification_RNASeq_DNASeq_Complete_2015_06_03.pdf

Sequencing summaries by lane and sample

	2015_05_28_RTSF_Output_DNA_S28-30_RNA_S1.3.5.7.9-11.13-15.17.19.22.23.26.xlsx

RTSF submission form

	RTSF_Submission_Form_RNASeq-Tumor-Controls_DNASeq-Controls_03-12-15.xlsx

Summary of flow cytometry measurements (Isolated CD4+ T-cell proportions per sample)

	2014_CD4_separation.pdf

Information about library used for sequencing (stranded)

	UG_Ovation_Single_Cell_RNA-Seq_System.pdf

Bioanalyzer measurements of genomic integrity from extracted RNA from each sample (some samples measured more than once)

	./bioanalyzer_results/

Explanation of samples, tissue, genomic material, and analyses

	Tumor_Sampling_Info.xlsx

