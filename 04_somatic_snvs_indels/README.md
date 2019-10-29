# Calling Somatic SNVs and Indels from Whole Genome Sequencing of Marek's Disease Lymphomas

summary

## Analyses in pipeline:

### Dependencies Used:
- MuSE (v1.0rc_c039ffa)
- MuTect (v1.1.7)
- JointSNVMix2 (v0.75)
- SomaticSniper (v1.0.5.0)
- VarDict (v1.4.4)
- VarScan2 (v2.4.1)
- Indelocator (v36.3336)
- LoFreq (2.1.2)
- Samtools (v0.1.12)
- bam-readcount (v0.7.4-unstable-39-dea4199)

### Formatting of Annotation and Reference:
- [The major contigs of Galgal5 and adjusted annotation](https://github.com/steepale/IKZF1_paper_code/tree/master/01_reference_prep)
- [The dbSNP set we converted to Galgal5](https://github.com/steepale/IKZF1_paper_code/tree/master/01_reference_prep/dbSNP_gg4-gg5-lift)

### Quality Control:

### Major Steps in Pipeline Analysis:
- Somatic SNVs were called with 6 algorithms
    1. MuSE
        - Default parameters
        - dbSNP file used to aid in filtering step
    2. MuTect
        - Default parameters
    3. JointSNVMix2 (EM Model)
        - Adjusted parameters:
            - Training performed with minimum coverage of 8 and 6 in normal and tumor samples, respectively.
        - Alternative alleles called if >= 0.95 probability
    4. SomaticSniper
        - SomaticSniper raw calls adjusted parameters:
            - Filtering reads with mapping quality less than 1
            - Filtering somatic snv output with somatic quality less than 20
            - Prior probability of a somatic mutation 0.01
        - Indel consensus sequences were called with Samtools:
            - Compute the consensus sequence
            - Print variants only
            - Only show lines/consensus with indels
        - SNPs were filtered with SomaticSniper (basic filter: default parameters)
        - SNPs underwent a second round of filtering at sites with:
            - Minimum base quality of 15
            - Minimum mapping quality of 1
    5. VarDict
        -  VarDict needs input bed files and [bed regions are recommended](https://github.com/AstraZeneca-NGS/VarDict/issues/2) to have 150 bp overlap for WGS data to call indels; bed files created with [create_bed_files_4_vardict.pl](https://github.com/hongenxu/MDV_proj/blob/master/somatic_snv_indel/create_bed_files_4_vardict.pl)
        - Variants were called with default parameters with one exception:
            - Variants with a variant allele frequency of 0.01 or greater were included
        - Variants filtered with a [python script](https://github.com/steepale/IKZF1_paper_code/blob/master/04_somatic_snvs_indels/vardict_fpfilter.py)
            - Total depth of position must be greater than or equal to 6
            - Variant depth must be greater than or equal to 2
            - Average mapping quality must be greater or equal to 20
            - The strand bias fisher p-value must be greater than or equal to 0.01
    6. VarScan2
        - Create named pipes (fifos) of samples
        - Collect Samtools mpileup format for VarScan2 somatic mutation calls:
            - Per-Base Alignment Quality disabled
            - Alignments with mapping quality score smaller than 1 were skipped
        - Somatic (SNPs and Indels) called with default parameters (germlines were also called)
        - Isolate Germline/LOH/Somatic calls from output with processSomatic (default parameters)

- Somatic Indels were called with 4 algorithms
    1. Indelocator
        - Indels called with default parameters
        - Filtered with [python script](https://github.com/steepale/IKZF1_paper_code/blob/master/04_somatic_snvs_indels/indelocator_fpfilter.py); Filtering parameters:
            - Total coverage greater than or equal to 6 at position in normal
            - Average mapping quality of consensus indel-supporting reads/reference-supporting reads greater than or equal to 20 in normal
            - Average quality of bases from consensus indel-supporting reads/from reference-supporting reads greater than or equal to 25 within NQS window in normal
            - Number of reads supporting consensus indel/any indel at the site greater than or equal to 2 in tumor
            - Total coverage greater than or equal to 6 at position in tumor
            - Average mapping quality of consensus indel-supporting reads/reference-supporting reads greater than or equal to 20 in tumor
            - Average quality of bases from consensus indel-supporting reads/from reference-supporting reads greater than or equal to 25 within NQS window in tumor
    2. VarDict
        - Described in SNV section
    3. VarScan2
        - Described in SNV section
    4. LoFreq
        - Indel qualities were added to bam files
        - Somatic indels were called with aid of dbSNP file and a minimum coverage of 6

## Data/Scripts:

## Genomic Datasets:

 
