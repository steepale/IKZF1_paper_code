#'---
#' title: "Examination of genotypes in Line 6, Line 7, and 6x7 F1s"
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
#' * Determine genotype similarity and differences between and among Line 6, Line 7, and F1 cross
#' * Generate a test and training dataset (SNPs only) to cross validate hard filting and mixed gaussian filtering

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
#BiocManager::install("R.utils")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr")
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

# Calculate the most common categorical value from vector or list
################################################################################
cat_mode <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
################################################################################

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Genotype calls (600K)
######################

# Load the 600K SNPs
K600_file <- "./data/truth_files/600K_Germline_SNPs_6x7_F1_Microarray.txt"
df <- read.table(file = K600_file, check.names = TRUE, header = TRUE, sep = '\t')
df <- as_tibble(df)

# Remove empty columns: X, X.1
df <- df %>% dplyr::select(-X, -X.1)
# Adjust the names
names(df) <- c("PROBE_SET_ID","CHR","POS","AFFY_SNP_ID","6_1","6_2","6_3","F1_1","F1_2","F1_3","F1_4","F1_5","F1_6","7_1","7_2","7_3")

# Clean up data
str(df)
summary(df)

#' ## Determine Genotype Similarirty amongst Line 6, Line 7, and F1 Genotypes
#'
#+ Geno Comp

################################################################################
#####     Determine Genotype Similarity      ###################################
################################################################################

# F1s
###########################

# Extract data frame of interest
#F1_df <- df %>% dplyr::select(PROBE_SET_ID,CHR,POS,AFFY_SNP_ID,F1_1, F1_2, F1_3, F1_4, F1_5, F1_6)

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan
F1_df <- df %>% dplyr::select(F1_1, F1_2, F1_3, F1_4, F1_5, F1_6)

# Convert rows into lists (20 seconds)
row_list <- F1_df %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
F1_df$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate sats about conserved and differing genotypes
F1_n_diff <- dim(F1_df %>% filter(GENO_FREQ < 1))[1]
F1_n_con <- dim(F1_df %>% filter(GENO_FREQ == 1))[1]
F1_n_all <- dim(F1_df)[1]
# Frequency of genotype conservation across 600K F1 genotypes
F1_freq_con = F1_n_con/F1_n_all

# 8885 sites differ out of 560086 in 6 samples
# 98.41% conserved

# F1s-3a (first 3 samples)
###########################

# Extract data frame of interest
#F1_df <- df %>% dplyr::select(PROBE_SET_ID,CHR,POS,AFFY_SNP_ID,F1_1, F1_2, F1_3, F1_4, F1_5, F1_6)

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan
F1_df_3a <- df %>% dplyr::select(F1_1, F1_2, F1_3)

# Convert rows into lists (20 seconds)
row_list <- F1_df_3a %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
F1_df_3a$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate sats about conserved and differing genotypes
F1_n_diff_3a <- dim(F1_df_3a %>% filter(GENO_FREQ < 1))[1]
F1_n_con_3a <- dim(F1_df_3a %>% filter(GENO_FREQ == 1))[1]
F1_n_all_3a <- dim(F1_df_3a)[1]
# Frequency of genotype conservation across 600K F1 genotypes
F1_freq_con_3a = F1_n_con_3a/F1_n_all_3a

# 6122 sites differ out of 560086 in 3 samples
# 98.90% conserved

# F1s-3b (last 3 samples)
###########################

# Extract data frame of interest
#F1_df <- df %>% dplyr::select(PROBE_SET_ID,CHR,POS,AFFY_SNP_ID,F1_1, F1_2, F1_3, F1_4, F1_5, F1_6)

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan
F1_df_3b <- df %>% dplyr::select(F1_4, F1_5, F1_6)

# Convert rows into lists (20 seconds)
row_list <- F1_df_3b %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
F1_df_3b$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate sats about conserved and differing genotypes
F1_n_diff_3b <- dim(F1_df_3b %>% filter(GENO_FREQ < 1))[1]
F1_n_con_3b <- dim(F1_df_3b %>% filter(GENO_FREQ == 1))[1]
F1_n_all_3b <- dim(F1_df_3b)[1]
# Frequency of genotype conservation across 600K F1 genotypes
F1_freq_con_3b = F1_n_con_3b/F1_n_all_3b

# 6996 sites differ out of 560086 in 3 samples
# 98.74% conserved

# Line 6
###########################

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan
L6_df <- df %>% dplyr::select('6_1', '6_2', '6_3')

# Convert rows into lists (20 seconds)
row_list <- L6_df %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
L6_df$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate the percentage of similarity
L6_df %>% filter(GENO_FREQ < 1)

# Calculate sats about conserved and differing genotypes
L6_n_diff <- dim(L6_df %>% filter(GENO_FREQ < 1))[1]
L6_n_con <- dim(L6_df %>% filter(GENO_FREQ == 1))[1]
L6_n_all <- dim(L6_df)[1]
# Frequency of genotype conservation across 600K F1 genotypes
L6_freq_con = L6_n_con/L6_n_all

# 748 sites differ out of 560086 in 3 samples
# 99.86% conserved

# Line 7
###########################

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan
L7_df <- df %>% dplyr::select('7_1', '7_2', '7_3')

# Convert rows into lists (20 seconds)
row_list <- L7_df %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
L7_df$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate the percentage of similarity
L7_df %>% filter(GENO_FREQ < 1)

# Calculate sats about conserved and differing genotypes
L7_n_diff <- dim(L7_df %>% filter(GENO_FREQ < 1))[1]
L7_n_con <- dim(L7_df %>% filter(GENO_FREQ == 1))[1]
L7_n_all <- dim(L7_df)[1]
# Frequency of genotype conservation across 600K F1 genotypes
L7_freq_con = L7_n_con/L7_n_all

# 4029 sites differ out of 560086 in 3 samples
# 99.27% conserved

# All samples: 3x Line 6, 3x Line 7, 6x F1s
###########################

# Generate rows for number and percent similarity between columns
# Note: Row-oriented work is difficult in R--https://resources.rstudio.com/webinars/thinking-inside-the-box-you-can-do-that-inside-a-data-frame-april-jenny-bryan
L67F1_df <- df %>% dplyr::select('6_1', '6_2', '6_3','7_1', '7_2', '7_3','F1_1','F1_2','F1_3','F1_4','F1_5','F1_6')

# Convert rows into lists (20 seconds)
row_list <- L67F1_df %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
L67F1_df$GENO_FREQ <- lapply(row_list,count_genos) %>% unlist()
# Calculate the percentage of similarity
L67F1_df %>% filter(GENO_FREQ < 1)

# Calculate sats about conserved and differing genotypes
L67F1_n_diff <- dim(L67F1_df %>% filter(GENO_FREQ < 1))[1]
L67F1_n_con <- dim(L67F1_df %>% filter(GENO_FREQ == 1))[1]
L67F1_n_all <- dim(L67F1_df)[1]
# Frequency of genotype conservation across 600K F1 genotypes
L67F1_freq_con = L67F1_n_con/L67F1_n_all

# 4029 sites differ out of 560086 in 3 samples
# 99.27% conserved

#' ## Combine Datasets and Save Results
#'
#+ Save Results

################################################################################
#####     Combine Datasets and Save Results      ###############################
################################################################################





# Bind GENO_FREQ columns into the original dataset and save
names(L7_df)[4] <- "L7_GENO_FREQ"
names(L6_df)[4] <- "L6_GENO_FREQ"
names(F1_df)[7] <- "F1_GENO_FREQ"

# Generate a tibble of data and conservation
df_all <- cbind(df,L6_df$L6_GENO_FREQ,L7_df$L7_GENO_FREQ,F1_df$F1_GENO_FREQ)
names(df_all)


