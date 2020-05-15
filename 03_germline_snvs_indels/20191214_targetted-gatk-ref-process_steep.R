#'---
#' title: "Process targetted GATK SNP calls from 20191214_gatk-600K-affy-calls_steep.sh"
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
#' * Process targetted GATK SNP calls from 20191214_gatk-600K-affy-calls_steep.sh for entry into 20191215_cross-val-germline_steep.R

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
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","ROCR","ff","GenomicRanges","BSgenome","rtracklayer", "caret")
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

# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * GATK 600K-Targetted Raw SNPs: unfiltered calls from 22 F1 control samples as g vcf (collective genotypes)
#' 
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# GATK 600K-Targetted Genotypes
######################

# Load the GATK 600K-targetted SNPs
df_file <- "./data/raw_snps_indels/collective-genotypes-600K-regions_steep.g.vcf.gz"
df_out <- str_replace(df_file, '.g.vcf.gz', '.g.vcf')
# Generate header
sys_cmd <- paste0('bgzip -d -c ',df_file, ' | grep "^#CHROM" | sed "s/#//" > ',df_out)
system(sys_cmd)
# Collect random data
sys_cmd <- paste0('bgzip -d -c ',df_file,' | grep -v "^#" >> ',df_out)
system(sys_cmd)

# Classes for columns
col_class <- c(rep('factor',1), 'integer', rep('factor',29))
# Read in file
GATK_K600_df <- read.table.ffdf(file = df_out, sep = '\t',header=TRUE,
                           colClasses=col_class, VERBOSE = TRUE, first.rows = 100000 ) %>% as_tibble()

# Examine data
dim(GATK_K600_df)
str(GATK_K600_df)
summary(GATK_K600_df)

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

# Save the sites to a file: will be used for targetted GATK analysis
write <- TRUE
if (write){
        print("file Saved!")
        K600_out <- K600_df %>% dplyr::select(CHROM,POS)
        K600_out$POS <- paste0(K600_out$POS,'-',K600_out$POS)
        # Write a list (for -L option in GATK)
        write.table(K600_out, './data/600K-regions_steep.list', sep =':',
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        # Write the entire dataframe
        write.table(K600_df, './data/600K-regions_steep.txt', sep ='\t',
                    quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Remove large memory items
rm(GG4to5_chain, K600_df_gr,K600_gr,K600_gr5,K600_df,K600_out)
