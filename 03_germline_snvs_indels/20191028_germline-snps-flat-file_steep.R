#'---
#' title: "Generation of a flat file for high-confidence Line 6 and Line 7 germline SNPs"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
#knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)



#' ## Explanation of Data 
#' Germline SNPs were called with GATK with this pipeline: [03_germline_snvs_indels](https://github.com/steepale/IKZF1_paper_code/tree/master/03_germline_snvs_indels).
#' 
#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Users/Alec/Documents/Bioinformatics/MDV_Project/germline_snps_indels/data/germline_snps/indv_samples/parents'
setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("R.utils")
#remove.packages("knitr")
#install.packages("knitr")

# Load dependencies
pacs...man <- c("dplyr","tibble","magrittr","stringr", "data.table","readr","R.utils","ggplot2","ggpubr")
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

################################################################################

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Line 6 and Line 7 SNPs from [germline SNP calling pipeline](https://github.com/steepale/IKZF1_paper_code/tree/master/03_germline_snvs_indels)
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Load the Line 6 SNPs
L6_file <- paste0(WD,"/002683_Line-6_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf")
L6 <- read.table(file = L6_file, comment.char = "", check.names = FALSE, header = TRUE, sep = '\t')
L6 <- as_tibble(L6)
# Adjust col name
colnames(L6)[1] <- 'CHROM'

# Load Line 7 data
L7_file <- paste0(WD,"/002684_Line-7_normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf")
L7 <- read.table(file = L7_file, comment.char = "", check.names = FALSE, header = TRUE, sep = '\t')
L7 <- as_tibble(L7)
# Adjust col name
colnames(L7)[1] <- 'CHROM'

#' ## Create flat file with the following format:
#' * Chromosome
#' * Position
#' * Line 6 Alternate Allele
#' * Line 7 Alternate Allele
#'
#+ Flat file

################################################################################
#####     Create Flat File      ################################################
################################################################################

# Select the columns of interest
L6 <- L6 %>% dplyr::select(CHROM, POS, REF, ALT)
L7 <- L7 %>% dplyr::select(CHROM, POS, REF, ALT)
# Rename the allese columns
colnames(L6)[3] <- "L6_REF"
colnames(L6)[4] <- "L6_ALLELE"
colnames(L7)[3] <- "L7_REF"
colnames(L7)[4] <- "L7_ALLELE"

# Perform a full-join
flat_data <- full_join(L6, L7, by = c("CHROM","POS"), copy = TRUE, suffix = c(".x", ".y"))
flat_data <- as_tibble(flat_data)

# Remove large objects
rm(L6)
rm(L7)

# Reformat for mutate function (factors are returned as ints)
flat_data$L6_REF <- as.character(flat_data$L6_REF)
flat_data$L6_ALLELE <- as.character(flat_data$L6_ALLELE)
flat_data$L7_REF <- as.character(flat_data$L7_REF)
flat_data$L7_ALLELE <- as.character(flat_data$L7_ALLELE)

# Replace NA values and add a column for reference
flat_data <- mutate(flat_data, L7_ALLELE = ifelse(is.na(L7_ALLELE), L6_REF, L7_ALLELE))
flat_data <- mutate(flat_data, L6_ALLELE = ifelse(is.na(L6_ALLELE), L7_REF, L6_ALLELE))
flat_data <- mutate(flat_data, REF = ifelse(is.na(L6_REF), L7_REF, L6_REF))
flat_data <- mutate(flat_data, REF = ifelse(is.na(L7_REF), L6_REF, L7_REF))

# Select the final columns
flat_data <- flat_data %>% dplyr::select(CHROM, POS, REF, L6_ALLELE, L7_ALLELE)

#' ## Visualize the output

# Visualize the output
str(flat_data)
head(flat_data)

# Save the flat file
flat_file <- paste0(WD,"/",date, "_L6-L7_SNPs-C99.9_",auth,".txt")
write.table(flat_data, file = flat_file, quote = FALSE, sep = '\t', row.names = FALSE)

#' ### Flat file saved as 20191028_L6-L7_SNPs-C99.9_steep.txt



