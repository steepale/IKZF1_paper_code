# Whole Genome Sequencing of Experimental Chicken F<sub>1</sub> Cross (6x7) from Lines 6<sub>3</sub> (Resistant to Marek's Disease) & 7<sub>2</sub> (Susceptible to Marek's Disease)

DNA was extracted from a blood sample of one line 6x7-F1 bird. All biological materials were prepared at the Avian Disease and Oncology Laboratory in East Lansing, Michigan. Whole genomic DNA from sample underwent WGS 150bp paired-end reads via the Illumina TruSeq Nano DNA Library Preparation Kit on Illumina HiSeq machines at the Michigan State University (MSU) RTSF Genomics Core (LIMS ID: CHE1458A17). Project designed and implemented by the Hans Cheng Lab (Hans.Cheng@ars.usda.gov).

## Supporting Annotation

Supporting annotation of birds, tumors, samples, and extracted DNA found in `supporting_documentation`

Summary of Submitted Samples to the MSU RTSF Genomics Core

	20131113_submitted-samples-msu-genomics_cheng.pdf

Sequencing summaries by lane and sample

	20140102_A_SeqProduction_Cheng.xlsx

## Sequencing Files

All raw FASTQ sequences are in `fastq_files`

	$ find fastq_files -name "*.fastq.gz"
	fastq_files/6x7-F1_GCCAAT_L001_R2_001.fastq.gz
    fastq_files/6x7-F1_GCCAAT_L002_R2_001.fastq.gz
    fastq_files/6x7-F1_GCCAAT_L001_R1_001.fastq.gz
    fastq_files/6x7-F1_GCCAAT_L002_R1_001.fastq.gz

## Data Integrity

MD5 file path: `fastq_files/line6x7-F1_blood.md5`

To perform MD5 check, from directory `fastq_files`

	$ md5sum -c line6x7-F1_blood.md5
	fastq_files/6x7-F1_GCCAAT_L001_R2_001.fastq.gz: OK
    fastq_files/6x7-F1_GCCAAT_L002_R2_001.fastq.gz: OK
    fastq_files/6x7-F1_GCCAAT_L001_R1_001.fastq.gz: OK
    fastq_files/6x7-F1_GCCAAT_L002_R1_001.fastq.gz: OK