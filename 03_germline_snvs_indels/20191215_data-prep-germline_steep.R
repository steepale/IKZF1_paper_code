#'---
#' title: "Predictive Model of Germline SNPs: Process Data"
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
#' * Prepare data for predictive modeling

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
BiocManager::install("BSgenome.Ggallus.UCSC.galGal4")
install.packages("earth")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","ROCR","ff","GenomicRanges","BSgenome","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal5","rtracklayer", "caret","pROC","modelr","ggplot2","e1071","doMC","glmnet", "lattice","MASS","pamr","pls","sparseLDA","lubridate","reshape2","kernlab", "klaR","latticeExtra","earth","partykit",
                "gtools")
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

# Function to load germline vcfs
################################################################################
load_normal_vcfs <- function(x, extract_gt = FALSE) {
        # Ensure 'x' is a string
        if ( class(x) != "character" ) {
                stop("'x' must be a string", class.= FALSE)
        }
        
        # Gunzip the file if necessary
        if (endsWith(x, '.gz')) {
                compressed <- TRUE
        }else{
                compressed <- FALSE 
        }
        
        if (compressed) {
                # Determine the number of lines to skip (with ##)
                system_cmd <- paste0('bgzip -d -c ',x,' | grep "^##" | wc -l')
                n_skip <- system(system_cmd, intern = TRUE) %>% as.numeric()
                # Determine the columns
                system_cmd <- paste0('bgzip -d -c ',x,' | grep "^#CHROM" | tr "\t" "\n" | wc -l')
                col_n <- system(system_cmd, intern = TRUE) %>% as.numeric()
        }else{
                # Determine the number of lines to skip (with ##)
                system_cmd <- paste0('grep "^##" ',x,'| wc -l')
                n_skip <- system(system_cmd, intern = TRUE) %>% as.numeric()
                # Determine the columns
                system_cmd <- paste0('grep "^#CHROM" ',x,'| tr "\t" "\n" | wc -l')
                col_n <- system(system_cmd, intern = TRUE) %>% as.numeric()
        }
        
        # Detemrine the column class for file loading
        col_class <- c("factor", "integer",rep("factor",3),"double",rep("factor", (col_n - 6) ))
        
        # Load vcf data
        #PL <- read.table(file = x, comment.char = '', skip = n_skip, check.names = FALSE, header = TRUE, sep = '\t') %>% as_tibble()
        PL <- read.table.ffdf(file = x, sep = '\t', skip=n_skip,header=TRUE,
                              colClasses=col_class, VERBOSE = TRUE, first.rows = 1000000 ) %>% as_tibble()
        
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        SAMPLE <- colnames(PL)[10]
        
        # Extract the genotypes if chosen
        if (extract_gt){
                # Extract the genotype information from the last column
                PL <- PL %>% separate(SAMPLE, c("GT", NA), sep = ':')
                # Select the columns of interest
                PL <- PL %>% dplyr::select(CHROM, POS, REF, ALT,GT)
                # Rename the allele columns
                colnames(PL)[5] <- paste0("GT")
                # Make GT factor
                PL$GT <- as.factor(PL$GT)
        }
        
        # Return the output
        PL
}
################################################################################

# Calculate the most common categorical value from vector or list
################################################################################
cat_mode_byrow <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
################################################################################

# Generate a frequency of genotypes
################################################################################
count_genos_byrow_F1 <- function(x) {
        # Ensure x is a list
        if (!is.list(x))
                stop("input is not a list")
        # Collect the most common genotype per row
        common <- paste0(x$GG5_REF,x$GG5_REF)
        # Count the common genotype occurance
        row_vector <- x %>% unlist()
        str_count(row_vector, common) %>% sum()
}
################################################################################

# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays with genotype similarities
#' * 600K Genotypes: training set for the VSQR model
#' * Line 6, Line 7, and 6x7 F1 germline variant (VSQR filtered) calls from WGS
#' * Line 6, Line 7, and 6x7 F1 germline variant (hard filtered) calls from WGS
#' * GATK 600K-Targetted Raw SNPs: unfiltered calls from 22 F1 control samples as g vcf (collective genotypes)
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
        train_file <- "./data/600K-regions-F1_steep.txt"
        train_out <- str_replace(train_file, '.txt', paste0('_',n,'.txt'))
        # Generate header
        sys_cmd <- paste0('head -n1 ',train_file, '> ',train_out)
        system(sys_cmd)
        # Collect random data
        sys_cmd <- paste0('grep -v "^#" ',train_file,' | shuf -n ',n,' >> ',train_out)
        system(sys_cmd)
}else{
        train_out <- "./data/600K-regions-F1_steep.txt"
}

# Classes for columns
col_class <- c('factor', 'integer', rep('factor',14), rep('double',7), rep('factor',7))
# Read in file
K600_df <- read.table.ffdf(file = train_out, sep = '\t',header=TRUE,
                           colClasses=col_class, VERBOSE = TRUE, first.rows = 100000 ) %>% as_tibble()

# Convert the Microarray ALT column from NA to '.'
K600_df <- K600_df %>%
        mutate(F1_ALT = ifelse(is.na(F1_ALT), '.', as.character(F1_ALT)))
K600_df$F1_ALT <- as.factor(K600_df$F1_ALT)

# Examine data
dim(K600_df)
str(K600_df)
summary(K600_df)

# Raw Germline SNP calls
######################

# Raw SNPs file
# Generate a file of randomly choosen rows for testing (consider putting 'if' statements in here if files already exist to save time)
n=99999
# Generate header
sys_cmd <- paste0('bgzip -d -c ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz | grep "^#" > ./data/raw_snps_indels/germline_raw_snps_indels_genotyped_',n,'.g.vcf')
system(sys_cmd)
# Collect random data
sys_cmd <- paste0('bgzip -d -c ./data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz | grep -v "^#" | shuf -n ',n,' >> ./data/raw_snps_indels/germline_raw_snps_indels_genotyped_',n,'.g.vcf')
system(sys_cmd)

# Load the raw SNP calls
raw_file <- './data/raw_snps_indels/germline_raw_snps_indels_genotyped.g.vcf.gz'
#raw_file <- paste0('./data/raw_snps_indels/germline_raw_snps_indels_genotyped_',n,'.g.vcf')
raw_df <- load_normal_vcfs(raw_file)

# Remove unneccessary columns from the start on attempt to reduce mem usage
raw_df <- raw_df %>% dplyr::select(-ID, -FILTER)

print("Seperating INFO Column")
# Transpose data frames and add annotation
raw_info <- raw_df %>% dplyr::select(INFO) %>% transpose()
# Search for ID, if there, generate column
#AC
INFO_COLUMNS <- c("AC","AF","AN","BaseQRankSum","ClippingRankSum","DP","ExcessHet","FS","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR")
# Iterate through each column and add the apprpriate annotation
for(x in INFO_COLUMNS){
        print(x)
        raw_df[[x]] <- str_match(raw_info, paste0(x,"=(.*?)(;|$)"))[,2] %>% as.double()
}

# Remove the INFO column
raw_df <- raw_df %>% dplyr::select(-INFO)

# Remove large memory items
rm(raw_info)

# Examine data
dim(raw_df)
str(raw_df)
summary(raw_df)

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

# Remove unneccessary columns from the start on attempt to reduce mem usage
GATK_K600_df <- GATK_K600_df %>% dplyr::select(-ID, -FILTER)

print("Seperating INFO Column")
# Transpose data frames and add annotation
GATK_K600_info <- GATK_K600_df %>% dplyr::select(INFO) %>% transpose()
# Search for ID, if there, generate column
#AC
INFO_COLUMNS <- c("AC","AF","AN","BaseQRankSum","ClippingRankSum","DP","ExcessHet","FS","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR")
# Iterate through each column and add the apprpriate annotation
for(x in INFO_COLUMNS){
        print(x)
        GATK_K600_df[[x]] <- str_match(GATK_K600_info, paste0(x,"=(.*?)(;|$)"))[,2] %>% as.double()
}

# Remove the INFO column
GATK_K600_df <- GATK_K600_df %>% dplyr::select(-INFO)

# Make QUAL numeric (DECISION)
GATK_K600_df <- GATK_K600_df %>%
        mutate(QUAL = ifelse(QUAL == '.', 0, as.double(QUAL)))

# Remove large memory items
rm(GATK_K600_info)

# Examine data
dim(GATK_K600_df)
str(GATK_K600_df)
summary(GATK_K600_df)

################################################################################
################# Combine Targetted and WGS SNP Calls ##########################
################################################################################

# T == GATK 600K-Targetted Genotypes (left)
# R == Raw Germline SNP calls (right)

dim(GATK_K600_df)
dim(raw_df)

# Get (R ∩ T) and R unique (right join)
df_rj <- GATK_K600_df %>% 
        dplyr::select(CHROM,POS,REF,ALT) %>% 
        right_join(raw_df, by = c('CHROM','POS','REF','ALT'))

# Remove most columns for sake of join
raw_laj <- raw_df %>% dplyr::select(CHROM,POS,REF,ALT)

# Get T unique (left anti join)
df_laj <- GATK_K600_df %>% 
        anti_join(raw_laj, by = c('CHROM','POS','REF','ALT'))

# Combine the dataframes into one dataframe
df_rt <- bind_rows(df_rj, df_laj) %>%
        arrange(CHROM,POS,REF,ALT)

# Exmaine dataset
dim(df_rt)
str(df_rt)
summary(df_rt)

# Remove large objects
rm(GATK_K600_df, raw_df, df_rj, df_laj)

################################################################################
########### Filter Genos to later make Training Data (600K SNP arrays) #########
################################################################################

# The frequencies of REF alleles were calculated to determine if microarray sites were called as SNPs (criteria: conserved across cohort)

# Remove any positions that: 
# - agree with reference 
# - have no SNP ID
# - are not conserved across cohort
train_tpfp <- K600_df %>% 
        dplyr::select('F1_1','F1_2','F1_3','F1_4','F1_5','F1_6','GG5_REF','AFFY_SNP_ID','F1_GENO_FREQ') %>% 
        filter(!is.na(AFFY_SNP_ID)) %>%
        filter(F1_GENO_FREQ == 1)
        
# Adjust genotype calls to characters
for( column in names(train_tpfp) ){
        train_tpfp[[column]] <- as.character(train_tpfp[[column]]) 
}

# Count the genotypes by row (4 minutes)
train_tpfp$REF_FREQ <- train_tpfp %>% f_pmap_aslist() %>% lapply(count_genos_byrow_F1) %>% unlist()
# Collect all calls in training set
train_all <- train_tpfp %>% dplyr::select(AFFY_SNP_ID, REF_FREQ)
# Determine whether the site is a SNP or not
train_all <- train_all %>% mutate(SNP_MA = ifelse(REF_FREQ <= 1, "1", "0"))

# Collect all calls in training set
#train_plot <- train_tpfp
# Determine whether the site is a SNP or not
#train_plot <- train_plot %>% mutate(SNP_MA = ifelse(REF_FREQ <= 1, "1", "0"))
# To make sure our filter is appropriate
#gp <- ggplot(data = train_plot) +
#        geom_jitter(
#                mapping = aes(x = REF_FREQ, y = F1_GENO_FREQ, color = SNP_MA, alpha = 0.001))
table(train_all %>% dplyr::select(REF_FREQ,SNP_MA))
# Subset again for sake of join
train_all <- train_all %>% dplyr::select(AFFY_SNP_ID, SNP_MA)

# Generate a revised train set
train_set_all <- inner_join(K600_df,train_all, by = "AFFY_SNP_ID")

# Generate a vector of homozygous alleles
homo_alleles <- c('AA','CC','TT','GG')
# Generate a vector of heterozygous alleles
prm <- gtools::permutations(n=4, r=2, v=c('A','T','C','G'))
hetero_alleles <- apply(prm, 1, function(x)paste0(x, collapse=''))

# Calculate F1_Genotype
train_set_all <- train_set_all %>%
        mutate(F1_GT = case_when(
                (F1_MODE %in% homo_alleles & F1_ALT == '.' & SNP_MA == '0') ~ '0/0',
                (F1_MODE %in% hetero_alleles & SNP_MA == '1') ~ '0/1',
                (F1_MODE %in% homo_alleles & F1_ALT != '.' & SNP_MA == '1')  ~ '1/1')
        )

# Data summary
dim(K600_df)
dim(train_set_all)

# Remove large items
rm(train_all)

################################################################################
################# Build Models of True Variant Prediction ######################
################################################################################

################################################################################
#### Dealing with Missing Values in Predictors (Cohort Level) ##################
################################################################################

# AN
##########################################
# A cohort count of allele number.
# Equals NA when all samples == './.' and 2 when atleast one samples demonstrates any genotype other than './.'
# GATK Source code supports calculation: GATKVariantContextUtils.java 
# Adjust AN NA values to represent 0 ()
df_rt <- df_rt %>% mutate(AN = ifelse(is.na(AN), 0, as.numeric(AN)))
df_rt$AN <- as.integer(df_rt$AN)

# AC
##########################################
# Allele count in genotypes, for each ALT allele, in the same order as listed
# Equals NA when there are no alt alleles (AKA only ./. or 0/0 across samples)
# GATK Source code supports calculation: GATKVariantContextUtils.java 

# Adjust AC NA values to represent 0 ()
df_rt <- df_rt %>% mutate(AC = ifelse(is.na(AC), 0, as.numeric(AC)))
df_rt$AC <- as.integer(df_rt$AC)

# AF
##########################################
# AF: Allele Frequency, for each ALT allele, in the same order as listed
# Equals NA when there are no alt alleles (AC == 0)
# GATK Source code supports calculation: GATKVariantContextUtils.java 
# Adjust AF NA values to represent 0 (occurs when AC == 0)
df_rt <- df_rt %>% mutate(AF = ifelse(is.na(AF) && AC == 0, 0, as.numeric(AF)))
df_rt$AF <- as.integer(df_rt$AF)

################################################################################
### Iterate through each sample and adjust dataset #############################
################################################################################

# Collect samples
samples <- names(df_rt)[startsWith(names(df_rt), 'X')]
# Generate empty dataframe for samples to be placed into
df_mod <- as_tibble()
# Model will be performed on a per-sample basis
for(sample in samples){
        
        ########################################################################
        ######################## Wrangle Data ##################################
        ########################################################################
        
        # Wrangle the FORMAT column
        #############################
        print(paste0('SAMPLE: ',sample))
        print('Wrangling Data...')
        # Other samples
        not_sample <- samples[samples %!in% sample]
        # Select columns needed
        sample_df <- df_rt %>% dplyr::select(-not_sample)
        # Separate the FORMAT and sample columns
        # Capture for format factors
        FRMT_FACTORS <- table(sample_df$FORMAT) %>% names()
        # Generate a tibble list
        tib_list <- list()
        # Iterate through each occurance of FORMAT column
        for( i in 1:length(FRMT_FACTORS) ){
                # Generate FACTOR
                FACTOR <- FRMT_FACTORS[i]
                # Generate the new columns
                new_cols <- str_split(FACTOR,':') %>% unlist()
                # Seperate the FACTOR column into new columns
                tib_list[[i]] <- sample_df %>% filter(FORMAT == FACTOR) %>% separate(sample, into = new_cols, sep = ':')
                # Make sure to add columns that do not exist
                all_names <- FRMT_FACTORS %>% str_split(':') %>% unlist() %>% unique()
                add_names <- all_names[all_names %!in% names(tib_list[[i]])]
                for (c in add_names) {
                        tib_list[[i]][[c]] <- NA
                }
                # Remove the FORMAT column
                tib_list[[i]] <- tib_list[[i]] %>% dplyr::select(-FORMAT)
        }
        # Combine the dataframes
        sample_df = do.call(rbind, tib_list)
        
        # Ensure columns are proper data types
        ############################
        
        # Seperate the AD column into REF_AD and ALT_AD
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        sample_df <- sample_df %>% separate(AD, into = c("REF_AD", "ALT_AD"), sep = ',')
        ## PL field will contain three numbers, corresponding to the three possible genotypes (0/0, 0/1, and 1/1). The PL values are "normalized" so that the PL of the most likely genotype (assigned in the GT field) is 0 in the Phred scale.
        sample_df <- sample_df %>% separate(PL, into = c("PL_00", "PL_01", "PL_11"), sep = ',')
        
        # Determine whether the site is a SNP or not (Totally Raw Calls)
        ############################
        sample_df <- sample_df %>% 
                mutate(SNP_GATK_RAW = ifelse(GT %in% c("./.","0/0"), '0', '1'))
        sample_df$SNP_GATK_RAW <- as.factor(sample_df$SNP_GATK_RAW)
        
        # Clean memory
        rm(tib_list)
        print('Finished Wrangling Data')
        
        ########################################################################
        ########## Subset GATK Call Set (Data for Training/Test Set) ###########
        ########################################################################
        # Collect only the variants that overlap with genotypes loci to create training set
        # Copy in training set
        train_df <- train_set_all %>% 
                dplyr::select(CHROM, POS, GG5_REF, F1_ALT, F1_MODE, F1_GT, SNP_MA)
        # Adjust column names
        colnames(train_df)[3] <- 'REF'
        #colnames(train_df)[4] <- 'ALT'
        
        # Perform a semi join to filter out variants from GATK that don't overlap with microarrays
        sample_df_gatk_ma <- inner_join(sample_df, train_df, by = c('CHROM','POS','REF'))
        #Some rows are called as SNP between WGS and Microarrays; however, they 
        #show different alternative alleles/genotypes: These data will be removed 
        # in training and testing sets (small n).
        rm_df <- sample_df_gatk_ma %>% 
                     filter(SNP_GATK_RAW == '1' & SNP_MA == '1' &
                                    (F1_ALT != ALT | F1_GT != GT))
        sample_df_gatk_ma <- anti_join(sample_df_gatk_ma, rm_df)
        
        # Examine the data
        #dim(sample_df_gatk_ma)
        #str(sample_df_gatk_ma)
        #summary(sample_df_gatk_ma)
        
        # Broadcast the sample
        sample_df_gatk_ma$SAMPLE <- as.factor(sample)
        
        # Add this dataframe to the final dataframe for model training/testing
        df_mod <- bind_rows(sample_df_gatk_ma, df_mod)
        
        # Clean up
        rm(train_df, rm_df)
}

########################################################################
######################## Pre-Model Evaluation ##########################
########################################################################

#' Prediction Accuracy (AKA success rate):
#' accuracy = (TP + TN)/(TP+TN+FP+FN)
#' 
#' The Error Rate (proportion of incorrectly classified examples):
#' error rate = (FP + FN)/(TP+TN+FP+FN) = 1 - accuracy
#' 
#' Kappa:
#' The kappa statistic adjusts accuracy by accounting for the possibility of a correct prediction by chance alone
#'Sensitivity: 
#'"The sensitivity of a model (... true positive rate), measures
#'the proportion of positive examples that were correctly classified."
#'
#'sensitivity = TP/(TP+FN)
#'
#'Specificity: 
#'"The specificity of a model (... true negative rate), measures
#'the proportion of negative examples that were correctly classified."
#'
#'specificity = TN/(TN+FP)
#'
#'Precision:
#'The precision (also known as the positive predictive value) is defines as the proportion of positive examples that are truly positive.
#'
#'precision = TP/(TP+FP)
#'
#'Recall:
#'A measure of how complete the results are. A model with high recall captures a large portion of the positive examples (wide breadth)
#'
#'recall = sensitiveity = TP/(TP+FN)
#'
#'~Lantz, B. Machine Learning with R. (Packt, 2019)

# Examine confusion matrix of unfiltered (EDA)
# Pos Pred Value == precision
df_mod$SNP_GATK_RAW <- as.factor(df_mod$SNP_GATK_RAW)
df_mod$SNP_MA <- as.factor(df_mod$SNP_MA)
confusionMatrix(df_mod$SNP_GATK_RAW,
                df_mod$SNP_MA, positive = "1")

#'The F-measure (F1-score, F-score):
#'A measure of model performance that combines precision and recall into a single number (harmonic mean).
#'
#'F-score = (2 x precision x recall) / (recall + precision)
#'
#'~Lantz, B. Machine Learning with R. (Packt, 2019)

# Precision and recall
prec <- posPredValue(df_mod$SNP_GATK_RAW,
                     df_mod$SNP_MA, positive = "1")
rec <- sensitivity(df_mod$SNP_GATK_RAW,
                   df_mod$SNP_MA, positive = "1")
# The F-score (harmonic mean)
f <- (2 * prec * rec) / (prec + rec)
f

#' Note: GATK raw calls are suprisingly robust

########################################################################
######################## Build Model of Data ###########################
########################################################################

########################################################################
################ Impute Variables at the Sample Level ##################
########################################################################

#' Features:
#' 
#' QUAL: phred-scaled quality score for the assertion made in ALT. i.e. -10log_10 prob(call in ALT is wrong). If ALT is ”.” (no variant) then this is -10log_10 p(variant), and if ALT is not ”.” this is -10log_10p(no variant). High QUAL scores indicate high confidence calls.
#' GT: Genotype
#' REF_AD: Allelic depths for the ref allele
#' ALT_AD: Allelic depths for the alt allele(s)
#' DP: Approximate read depth (reads with MQ=255 or with bad mates are filtered)
#' GQ: Genotype Quality
#' PGT: Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another
#' PID: Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group
#' PL_00: Phred-scaled likelihoods for (00) genotype
#' PL_01: Phred-scaled likelihoods for (01) genotype
#' PL_11: Phred-scaled likelihoods for (11) genotype
#' AC: Allele count in genotypes, for each ALT allele, in the same order as listed
#' AF: Allele Frequency, for each ALT allele, in the same order as listed
#' AN: Total number of alleles in called genotypes
#' BaseQRankSum: Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities
#' ClippingRankSum: Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases
#' ExcessHet: Phred-scaled p-value for exact test of excess heterozygosity
#' FS: Phred-scaled p-value using Fisher's exact test to detect strand bias
#' InbreedingCoeff: Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation
#' MLEAC: Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed
#' MLEAF: Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed
#' MQ: RMS Mapping Quality
#' MQRankSum: Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities
#' QD: Variant Confidence/Quality by Depth
#' ReadPosRankSum: Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias
#' SOR: Symmetric Odds Ratio of 2x2 contingency table to detect strand bias
#' RGQ: Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)

# Ensure columns are proper data types
############################

# Iterate through each column of tibble and if character, adjust to factor
df_mod <- df_mod %>% mutate_if(sapply(df_mod, is.character), as.factor)
# Determine columns that should be integers and convert them
c2int <- c('REF_AD', 'ALT_AD','DP','GQ','RGQ','PL_00', "PL_01", "PL_11", "AC",'MLEAC')
for (c2i in c2int) {
        df_mod[[c2i]] <- as.integer(df_mod[[c2i]])
}
# Determine columns that should be doubles and convert them
c2dbl <- c('AF','AN','ExcessHet','FS','InbreedingCoeff','MLEAF','MQ','QD','SOR')
for (c2d in c2dbl) {
        df_mod[[c2d]] <- as.numeric(df_mod[[c2d]])
}
        
# Examine Data
str(df_mod)
summary(df_mod)

# Write the dataframe to file
out_file <- './data/processed-data_steep.txt'
write.table(df_mod, file = out_file, quote = FALSE, sep = '\t')

sessionInfo()
