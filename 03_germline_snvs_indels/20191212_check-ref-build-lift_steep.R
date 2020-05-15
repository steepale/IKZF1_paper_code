#'---
#' title: "Check Reference for Datasets and LiftOver when Necessary"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels'
setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
# BiocManager::install("ade4")
#install.packages("ade4")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","readxl", "GenomicRanges","BSgenome","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal5","rtracklayer")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

############################################################
##### Functions ############################################
############################################################

# Make the 'not in' operator
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}

# Capture the Date
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "steep"

## explicit gc, then execute `expr` `n` times w/o explicit gc, return timings (~"Jenny Bryan, updating work of Winston Chang")
benchmark <- function(n = 1, expr, envir = parent.frame()) {
        expr <- substitute(expr)
        gc()
        map(seq_len(n), ~ system.time(eval(expr, envir), gcFirst = FALSE))
}

# Generate a frequency of genotypes
################################################################################
count_genos <- function(x) {
        # Ensure x is a list
        if (!is.list(x))
                stop("input is not a list")
        # Collect the most common genotype per row
        common <- cat_mode(x) %>% unlist() %>% as.character()
        # Count the common genotype occurance
        row_vector <- x %>% unlist()
        sum(str_count(row_vector, common)) / (length(row_vector)-sum(str_count(row_vector,'---')))
}
################################################################################

# Calculate the most common categorical value from vector or list
################################################################################
cat_mode <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
################################################################################

# Function to load germline vcfs
################################################################################
load_normal_vcfs <- function(x) {
        # Ensure 'x' is a string
        if ( class(x) != "character" ) {
                stop("'x' must be a string", class.= FALSE)
        }
        # Assign the vcf file
        ND <- "/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/germline_snps/indv_samples/parents"
        vcf.gz <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf.gz')
        vcf <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf')
        # Gunzip the file
        gunzip(vcf.gz)
        # Load vcf data
        PL <- read.table(file = vcf, comment.char = "", check.names = FALSE, header = TRUE, sep = '\t')
        PL <- as_tibble(PL)
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        SAMPLE <- colnames(PL)[10]
        # gzip the file
        gzip(vcf)
        
        # Extract the genotype information from the last column
        PL <- PL %>% separate(SAMPLE, c("GT", NA), sep = ':')
        # Select the columns of interest
        PL <- PL %>% dplyr::select(CHROM, POS, REF, ALT,GT)
        
        # Rename the allele columns
        colnames(PL)[5] <- paste0("GT")
        
        # Return the output
        as_tibble(PL)
}
################################################################################

# Function to load germline vcfs (total)
################################################################################
load_total_vcf <- function(x) {
        # Ensure 'x' is a string
        if ( class(x) != "character" ) {
                stop("'x' must be a string", class.= FALSE)
        }
        # Assign the vcf file
        ND <- "/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/germline_snps/indv_samples/parents"
        vcf.gz <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf.gz')
        vcf <- paste0(ND,'/',x,'_normals_recal_famprioirs_gfiltered_denovo_germline_99.9_10K.g.vcf')
        # Gunzip the file
        gunzip(vcf.gz)
        # Load vcf data
        PL <- read.table(file = vcf, comment.char = "", check.names = FALSE, header = TRUE, sep = '\t')
        PL <- as_tibble(PL)
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        SAMPLE <- colnames(PL)[10]
        # gzip the file
        gzip(vcf)
        # Return the output
        as_tibble(PL)
}
################################################################################

# Function to load vcfs
################################################################################
load_vcfs <- function(x) {
        # Ensure 'x' is a string
        if ( class(x) != "character" ) {
                stop("'x' must be a string", class.= FALSE)
        }
        # Gunzip the file if necessary
        if (endsWith(x, '.gz')) {
                gunzip(x)
                x <- str_replace(x, '.gz', '')
        }
        
        # Determine the number of lines to skip (with ##)
        system_cmd <- paste0('grep "^##" ',x,' | wc -l')
        n_skip <- system(system_cmd, intern = TRUE) %>% as.numeric()
        
        # Load vcf data
        PL <- read.table(file = x, comment.char = '', skip = n_skip, check.names = FALSE, header = TRUE, sep = '\t') %>% as_tibble(PL)
        
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        
        # Determine if file needs to be gzipped
        if (endsWith(x, '.vcf')) {
                gzip(x)
                x <- str_replace(x, '.vcf', '.vcf.gz')
        }
        
        # Return the output
        PL
}
################################################################################

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays
#' 
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Chain file for liftover
######################
# Install a chain file for liftover
# Chain format: https://genome.ucsc.edu/goldenpath/help/chain.html
# Liftover file downloaded (20191205) from: http://hgdownload.soe.ucsc.edu/goldenPath/galGal4/liftOver/

# Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays with genotype similarities
######################

# Load the 600K SNPs
K600_file <- "./data/20191212_L67F1-genotype-frequencies_steep.txt"
K600_df <- read.table(file = K600_file, check.names = FALSE, header = TRUE, sep = '\t') %>% as_tibble()

# Examine data
str(K600_df)
summary(K600_df)

# 600K Genotypes: training set for the VSQR model
######################

# Load the 600K training SNPs
train_file <- "./data/truth_files/collective_samples_600K_shared_snps_formated.g.vcf.gz"
train_df <- load_vcfs(train_file)
# Select columns
train_df <- train_df %>% dplyr::select(CHROM,POS,REF,ALT)

# Examine data
str(train_df)
summary(train_df)

#' ## Determine the reference build
#'
#+ Determine reference

################################################################################
#####     Determine the Reference Build      ###################################
################################################################################

# Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays with genotype similarities
##############

# 600K SNP calls previously used to train the VSQR model
# Select columns of interest from K600
K600_F1 <- K600_df %>% filter(F1_GENO_FREQ == 1) %>% dplyr::select(CHR,POS,F1_1,F1_GENO_FREQ)

# Adjust the col name of K600
names(K600_F1)[1] <- 'CHROM'

# Remove genotypes with no chromosome
K600_F1 <- K600_F1 %>% filter(CHROM != '---')

# Take a random sample of df
K600_F1_10K <- sample_n(K600_F1, 10000)

# Create a start and end column because getSeq might otherwise error
K600_F1_10K$END <- K600_F1_10K$POS

# Obtain the reference sequences (getSeq is vectorized)
# Galgal4
K600_F1_10K$GG4_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal4::Ggallus, names = paste0('chr',K600_F1_10K$CHROM), start=K600_F1_10K$POS, end=K600_F1_10K$POS,
                         strand="+", as.character=TRUE)
# Galgal5
K600_F1_10K$GG5_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal5::Ggallus, names = paste0('chr',K600_F1_10K$CHROM), start=K600_F1_10K$POS, end=K600_F1_10K$POS,
                         strand="+", as.character=TRUE)

K600_F1_10K <- K600_F1_10K %>% 
        filter(POS %!in% c(11153298, 11096832,11080108,14262517,14294941,6286730)) %>% 
        filter(CHROM != 'LGE22C19W28_E50C23')

# Make sure all reference positions are true 
row_list <- K600_F1_10K %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
F1_df_3a$GENO_FREQ <- lapply(row_list,vectorized_grepl) %>% unlist()
eow_list <- row_list[[1]]
# Generate a frequency of genotypes
################################################################################
vectorized_grepl <- function(row_list) {
        # Ensure row_list is a list
        if (!is.list(row_list))
                stop("input is not a list")
        # determine if they match
        
}
################################################################################

table(K600_F1_10K$GG4_REF)

table(grepl(pattern = K600_F1_10K$GG4_REF[1],K600_F1_10K$F1_1))
table(grepl(pattern = K600_F1_10K$GG5_REF[1],K600_F1_10K$F1_1))

# Test GG4
K600_F1_10K[str_count(K600_F1_10K$F1_1, pattern = K600_F1_10K$GG4_REF) >= 1,]
K600_F1_10K[str_count(K600_F1_10K$F1_1, pattern = K600_F1_10K$GG4_REF) == 0,]
# Test GG5
K600_F1_10K[str_count(K600_F1_10K$F1_1, pattern = K600_F1_10K$GG5_REF) >= 1,]
K600_F1_10K[str_count(K600_F1_10K$F1_1, pattern = K600_F1_10K$GG5_REF) == 0,]

# Conclusion: Reference build is GG4 and all REFs refer to "+" strand
