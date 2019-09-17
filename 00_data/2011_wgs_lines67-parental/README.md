# Whole Genome Sequencing of Experimental Chicken Lines 6<sub>3</sub> (Resistant to Marek's Disease) & 7<sub>2</sub> (Susceptible to Marek's Disease)

DNA was extracted from blood samples of 6 male birds from each line and then pooled by line--DNA from six males of line 6 and six males of line 7 make the line 6 and line 7 samples, respectively. All biological materials were prepared at the Avian Disease and Oncology Laboratory in East Lansing, Michigan. Sequencing was performed in 2011 by DNA Landmarks (Montreal, Canada) using an Illumina HiSeq (100bp paired-end reads).

## Sequencing Files

All raw FASTQ sequences are in `fastq_files/illumina`

	$ find fastq_files/illumina -name "*.fastq.gz"
	fastq_files/illumina/002683_Line-6_ACTTGA_L001_R2_001.fastq.gz
    fastq_files/illumina/002683_Line-6_ACTTGA_L002_R1_001.fastq.gz
    fastq_files/illumina/002683_Line-6_ACTTGA_L001_R1_001.fastq.gz
    fastq_files/illumina/002684_Line-7_GATCAG_L003_R1_001.fastq.gz
    fastq_files/illumina/002683_Line-6_ACTTGA_L002_R2_001.fastq.gz
    fastq_files/illumina/002684_Line-7_GATCAG_L001_R1_001.fastq.gz
    fastq_files/illumina/002684_Line-7_GATCAG_L003_R2_001.fastq.gz
    fastq_files/illumina/002684_Line-7_GATCAG_L002_R2_001.fastq.gz
    fastq_files/illumina/002683_Line-6_ACTTGA_L003_R2_001.fastq.gz
    fastq_files/illumina/002684_Line-7_GATCAG_L002_R1_001.fastq.gz
    fastq_files/illumina/002683_Line-6_ACTTGA_L003_R1_001.fastq.gz
    fastq_files/illumina/002684_Line-7_GATCAG_L001_R2_001.fastq.gz

## Data Integrity

MD5 file path: `fastq_files/illumina/line6_line7_illumina.md5`

To perform MD5 check, from directory `fastq_files/illumina`

	$ md5sum -c line6_line7_illumina.md5
	002683_Line-6_ACTTGA_L001_R2_001.fastq.gz: OK
    002683_Line-6_ACTTGA_L002_R1_001.fastq.gz: OK
    002683_Line-6_ACTTGA_L001_R1_001.fastq.gz: OK
    002684_Line-7_GATCAG_L003_R1_001.fastq.gz: OK