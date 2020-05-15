#'---
#' title: "LiftOver 600K sites from GG4 to GG5"
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

#' # Goals of Analysis
#' * LiftOver 600K sites from GG4 to GG5

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
#BiocManager::install("ff")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","ROCR","ff","GenomicRanges","BSgenome","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal5","rtracklayer", "caret")
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

## explicit gc, then execute `expr` `n` times w/o explicit gc, return timings
benchmark <- function(n = 1, expr, envir = parent.frame()) {
        expr <- substitute(expr)
        gc()
        map(seq_len(n), ~ system.time(eval(expr, envir), gcFirst = FALSE))
}

# Calculate the most common categorical value from vector or list
################################################################################
cat_mode <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
################################################################################

# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################

# Determine the most common allele amongst genotypes (Similar to "cat_mode")
################################################################################
alt_mode <- function(x, ref_col, geno_col) {
        # Ensure x is a list
        if (!is.list(x))
                stop("input is not a list")
        
        # Collect geno alleles (split)
        alleles <- x[[geno_col]] %>% str_split(pattern = '') %>% unlist()
        
        # If condition
        bool_con <- all(unique(alleles) %in% x[[ref_col]])
        
        # Collect the alternative allele, if there is one, that isn't in the ref allele
        if (bool_con){
                ALT <- NA
        }else{
                ALT <- alleles[alleles %!in% x[[ref_col]]] %>% unique()
                ## If more than 2 alleles, condense into 1 string
                if (length(ALT) > 1){
                        ALT <- paste(ALT, collapse = ',')
                }
        }
        # The output
        ALT
}
################################################################################


#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * 600K Genotypes: training set for the machine learning model
#' 
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# 600K Genotypes
######################
# Decision to perform random sub sampling
subset <- FALSE
if(subset) {
        # Generate a file of randomly choosen values for testing
        n=99999
        # Load the 600K training SNPs
        train_file <- "./data/20191212_L67F1-genotype-frequencies_steep.txt"
        train_out <- str_replace(train_file, '.txt', paste0('_',n,'.txt'))
        # Generate header
        sys_cmd <- paste0('head -n1 ',train_file, '> ',train_out)
        system(sys_cmd)
        # Collect random data
        sys_cmd <- paste0('grep -v "^#" ',train_file,' | shuf -n ',n,' >> ',train_out)
        system(sys_cmd)
}else{
        train_out <- "./data/20191212_L67F1-genotype-frequencies_steep.txt"
}

# Classes for columns
col_class <- c(rep('factor',2), 'integer', rep('factor',13), rep('double',7))
# Read in file
K600_df <- read.table.ffdf(file = train_out, sep = '\t',header=TRUE,
                           colClasses=col_class, VERBOSE = TRUE, first.rows = 100000 ) %>% as_tibble()

# Adjust the Column names
names(K600_df)[2] <- 'CHROM'

# Examine data
dim(K600_df)
str(K600_df)
summary(K600_df)

#' ## Perform liftover from Galgal4 to Galgal5 on 600K Affy Loci
#'
#+ GG4 to GG5

################################################################################
#####     Perform Liftover      ################################################
################################################################################

# Adjust the name of Chromosomes Z and X before entering into GRanges object
K600_df_gr <- K600_df
K600_df_gr$CHROM <- K600_df$CHROM %>% as.character()
K600_df_gr[K600_df_gr$CHROM == 'Z',]$CHROM <- 'chrZ'
K600_df_gr[K600_df_gr$CHROM == 'LGE64',]$CHROM <- 'chrLGE64'
K600_df_gr$CHROM <- K600_df_gr$CHROM %>% as.factor()

# Remove unneccessary CHROM
K600_df_gr <- K600_df_gr %>% filter(CHROM != "LGE22C19W28_E50C23") %>% filter(CHROM != "---")
K600_df_gr$CHROM <- as.character(K600_df_gr$CHROM) %>% as.factor()

# Load the data into a GRanges object
K600_gr <- makeGRangesFromDataFrame(K600_df_gr,
                                    keep.extra.columns = TRUE,
                                    ignore.strand=FALSE,
                                    seqnames.field='CHROM',
                                    start.field = 'POS',
                                    end.field = 'POS',
                                    strand.field = "+")
# Add 'chr' to CHROM names
seqlevelsStyle(K600_gr) = "UCSC"
seqlevels(K600_gr)
#seqnames(K600_gr[seqnames(K600_gr) == 'Z']) <- 'chrZ'

# Decompress the chain file
gunzip('./data/galGal4ToGalGal5.over.chain.gz')
# Import the Chain file
GG4to5_chain <- import.chain('./data/galGal4ToGalGal5.over.chain')
# Compress the chain file
gzip('./data/galGal4ToGalGal5.over.chain')

# Perform liftover
K600_gr5 <- liftOver(K600_gr, GG4to5_chain)
K600_gr5 <- unlist(K600_gr5)
# Convert to tibble
K600_df <- as_tibble(K600_gr5)
# Drop old REF columns
K600_df <- K600_df %>% dplyr::select(-end, -strand, -width)

# Adjust names
names(K600_df)[1:2] <- c('CHROM','POS')

# Collect the GG5 References
K600_df$GG5_REF <- getSeq(BSgenome.Ggallus.UCSC.galGal5::Ggallus, names = K600_df$CHROM, start=K600_df$POS, end=K600_df$POS,strand="+", as.character=TRUE)

# Remove the 'chr' from CHROM field
K600_df$CHROM <- str_replace_all(K600_df$CHROM, 'chr', '')

# Dimensions
dim(K600_df)
#K600_df <- K600_df2
#K600_df <- head(K600_df, n = 20000)
#x <- K600_df[824,]

#K600_df %>% dplyr::select(X6_MODE, GG5_REF, X6_ALT) %>% View()
# Add 3 columns of Annotation: L6_MODE, F1_MODE, L7_MODE
# Collect the appropriate ALT Approximation per row (11 minutes total)
K600_df$F1_MODE <- K600_df %>% 
        dplyr::select(F1_1,F1_2,F1_3,F1_4,F1_5,F1_6) %>%
        f_pmap_aslist() %>%
        lapply(cat_mode) %>% unlist()
K600_df$X6_MODE <- K600_df %>% 
        dplyr::select(X6_1, X6_2, X6_3) %>%
        f_pmap_aslist() %>%
        lapply(cat_mode) %>% unlist()
K600_df$X7_MODE <- K600_df %>% 
        dplyr::select(X7_1, X7_2, X7_3) %>%
        f_pmap_aslist() %>%
        lapply(cat_mode) %>% unlist()
# Collect the alternative allele (if any) from the most common genotype
K600_df$F1_ALT <- K600_df %>% 
        dplyr::select(F1_MODE, GG5_REF) %>%
        f_pmap_aslist() %>%
        lapply(alt_mode, ref_col = "GG5_REF", geno_col = "F1_MODE") %>% unlist() %>% as.factor()
K600_df$X6_ALT <- K600_df %>% 
        dplyr::select(X6_MODE, GG5_REF) %>%
        f_pmap_aslist() %>%
        lapply(alt_mode, ref_col = "GG5_REF", geno_col = "X6_MODE") %>% unlist() %>% as.factor()
K600_df$X7_ALT <- K600_df %>% 
        dplyr::select(X7_MODE, GG5_REF) %>%
        f_pmap_aslist() %>%
        lapply(alt_mode, ref_col = "GG5_REF", geno_col = "X7_MODE") %>% unlist() %>% as.factor()

## Select columns and filter out variants not conserved across all F1's
K600_df_F1 <- K600_df %>% filter(F1_GENO_FREQ == 1)
## Select columns and filter out variants not conserved across all F1's
K600_df_X6 <- K600_df %>% filter(L6_GENO_FREQ == 1)
## Select columns and filter out variants not conserved across all F1's
K600_df_X7 <- K600_df %>% filter(L7_GENO_FREQ == 1)

# Save the sites to a file: will be used for targetted GATK analysis
# F1:
###################
K600_out <- K600_df_F1 %>% dplyr::select(CHROM,POS)
K600_out$POS <- paste0(K600_out$POS,'-',K600_out$POS)
# Write a list (for -L option in GATK)
write.table(K600_out, './data/600K-regions-F1_steep.list', sep =':',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# Write the entire dataframe
write.table(K600_df, './data/600K-regions-F1_steep.txt', sep ='\t',
            quote = FALSE, row.names = FALSE, col.names = TRUE)
# Line 6:
###################
K600_out <- K600_df_X6 %>% dplyr::select(CHROM,POS)
K600_out$POS <- paste0(K600_out$POS,'-',K600_out$POS)
# Write a list (for -L option in GATK)
write.table(K600_out, './data/600K-regions-L6_steep.list', sep =':',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# Write the entire dataframe
write.table(K600_df, './data/600K-regions-L6_steep.txt', sep ='\t',
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Line 7:
###################
K600_out <- K600_df_X7 %>% dplyr::select(CHROM,POS)
K600_out$POS <- paste0(K600_out$POS,'-',K600_out$POS)
# Write a list (for -L option in GATK)
write.table(K600_out, './data/600K-regions-L7_steep.list', sep =':',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# Write the entire dataframe
write.table(K600_df, './data/600K-regions-L7_steep.txt', sep ='\t',
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Remove large memory items
rm(GG4to5_chain, K600_df_gr,K600_gr,K600_gr5,K600_df,K600_out)
