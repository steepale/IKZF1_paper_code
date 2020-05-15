#'---
#' title: "K-Folds Cross Validation of Germline SNP Calling (VSQR & Hard Filters)"
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
#' * Perform k fold cross validation on germline snp calls (hard filtered and mixed gaussian)

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
                              colClasses=col_class, VERBOSE = TRUE, first.rows = 100000 ) %>% as_tibble()
        
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        SAMPLE <- colnames(PL)[10]
        
        # Extract the genotypes if choosen
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
raw_file <- paste0('./data/raw_snps_indels/germline_raw_snps_indels_genotyped_',n,'.g.vcf')
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

purities_df
ggplot

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

#summary(df_mod$AC)
#is.na(df_mod$AC) %>% table()
#df_mod %>% filter(is.na(AC)) %>% 
#        filter(GT %!in% c('./.', '0/0')) %>% 
#        dplyr::select(SNP_MA)
#View(df_mod %>% filter(is.na(AC)))
#View(df_mod %>% filter(AC < 2))

# Custom adjustment of variables (Sample-Context Dependent)
######################################################

# ALT_AD (Per sample)
##############################

# Is a value of NA for ALT_AD indicative of no allele detected? -- YES
df_mod %>% filter(is.na(ALT_AD)) %>% dplyr::select(ALT) %>% table()
df_mod %>% filter(is.na(ALT_AD)) %>% dplyr::select(GT) %>% table()
# Adjust ALT_AD NA values to represent 0
df_mod <- df_mod %>% mutate(ALT_AD = ifelse(is.na(ALT_AD), 0, as.numeric(ALT_AD)))
df_mod$ALT_AD <- as.integer(df_mod$ALT_AD)

# Consider Changing GQ to 0 if NA
#df_mod <- df_mod %>% mutate(GQ = ifelse(is.na(GQ), 0, as.numeric(GQ)))
#df_mod$GQ <- as.numeric(df_mod$GQ)


# Investigation of varibales with models and plots
######################################################

#df_mod %>% filter(AC < 1) %>% View()
summary(df_mod$AC)
table(df_mod$AN)

"""

# Boxplot of AN vs GT
ggplot(data = df_mod, 
       mapping = aes(x = GT, y = ALT_AD, colour = SNP_GATK_RAW)) +
        geom_boxplot() +
        geom_jitter(height = 0, alpha = 0.4)
        #scale_y_continuous(breaks = 1:10) +
        #ylim(0,10)
        
        df_mod %>% filter(GT == '0/0') %>% 
                dplyr::select(AC) %>% 
                unlist() %>% 
                hist(breaks = 40)
        
        sim1_mod <- lm(REF_AD ~ GQ, data = df_mod)
        summary(sim1_mod)
        grid <- df_mod %>% data_grid(GQ)
        grid <- grid %>% add_predictions(sim1_mod)
        ggplot(df_mod, aes(GQ)) +
                geom_point(
                        aes(y = REF_AD,
                            colour = GT),
                        alpha = 0.2,
                        ) +
                geom_line(
                        aes(y = pred),
                        data = grid,
                        colour = 'red',
                        size = 1
                )
        
        df_mod %>% filter(GQ >= 80) %>% 
                dplyr::select(QUAL) %>% 
                unlist() %>% 
                hist()
        
        sim1_mod <- lm(QUAL ~ GQ, data = df_mod)
        summary(sim1_mod)
        grid <- df_mod %>% data_grid(GQ)
        grid <- grid %>% add_predictions(sim1_mod)
        ggplot(df_mod, aes(GQ)) +
                geom_point(
                        aes(y = QUAL),
                        alpha = 0.08) +
                geom_line(
                        aes(y = pred),
                        data = grid,
                        colour = 'red',
                        size = 1
                )
        # Addition plots of model
        plot(sim1_mod)
        
        
        #relabel certain columns to ordinal
        #ordinals <- c('NUM_TOOLS','FILTER')
        #for (ord in ordinals) {
        #        df[[ord]] <- ordered(as.factor(df[[ord]])) 
        #}
"""

# Examine/Collect Variables
###############################
str(df_mod)
summary(df_mod)

########################################################################
#################### Hard-Filter Evaluation ############################
########################################################################

# Hard Filter Parameters
# QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0
        
        df_mod %>%
                filter(!is.na(QD)) %>%
                dplyr::select(SNP_MA) %>%
                table()
        
        df_mod %>%
                filter(QD > 2) %>%
                filter(FS < 60) %>%
                filter(MQ > 40) %>%
                filter(MQRankSum > -12.5) %>%
                filter(ReadPosRankSum > -8) %>%
                dplyr::select(SNP_MA) %>%
                table()
        
        # Dependents:
        # 'QUAL': Missed calls in 0 - 0.5
        # 'GT'
        # 'REF_AD'
        # 'ALT_AD': Missed calls in 0 - 0.5
        # 'DP'
        # 'AN'
        # 'AF'
        
        # Perform Hard Filtering
        summary(df_mod$AF)
        ggplot(df_mod, aes(AF, fill = SNP_MA, colour = SNP_MA)) +
                geom_density(alpha = 0.2, na.rm = TRUE) +
                xlim(0.6,1.4)
        
        df_mod <- df_mod %>%
                mutate(SNP_GATK_HARD = ifelse( (QUAL >= 0.5) , '1', '0'))
        df_mod$SNP_GATK_HARD <- as.factor(df_mod$SNP_GATK_HARD)
        
        # Examine the confuson matrix
        confusionMatrix(df_mod$SNP_GATK_RAW,
                        df_mod$SNP_MA, positive = "1")
        confusionMatrix(df_mod$SNP_GATK_HARD,
                        df_mod$SNP_MA, positive = "1")
        
        #'The F-measure (F1-score, F-score):
        #'A measure of model performance that combines precision and recall into a single number (harmonic mean).
        #'
        #'F-score = (2 x precision x recall) / (recall + precision)
        #'
        #'~Lantz, B. Machine Learning with R. (Packt, 2019)
        
        # Precision and recall
        prec <- posPredValue(sample_df_gatk_ma$SNP_GATK_RAW,
                             sample_df_gatk_ma$SNP_MA, positive = "1")
        rec <- sensitivity(sample_df_gatk_ma$SNP_GATK_RAW,
                           sample_df_gatk_ma$SNP_MA, positive = "1")
        # The F-score (harmonic mean)
        f <- (2 * prec * rec) / (prec + rec)
        f
        
        #' Note: GATK raw calls are suprisingly robust
        
        # Start with the variables intended to be used in the model
        drop_col <- c('CHROM','POS','REF','ALT','AFFY_SNP_ID','VAL','PID')
        dependents <- names(df_mod)[names(df_mod) %!in% drop_col]
        
        
        # Correct for NA values (Categorical Variables)
        #############################
        # Drop dependents with too many NAs
        # Create a function to add to drop columns if 20% or more of variables are NA
        na_num <- c()
        na_cat <- c()
        for (var in dependents){
                # Determine the amount of NA
                n_na <- table(is.na(df[[var]]))[2] %>% as.numeric()
                # Condition to turn n_na to 0 if it equals NA
                if (is.na(n_na)) n_na <- 0
                n_not_na <- table(is.na(df[[var]]))[1] %>% as.numeric()
                freq_na <- n_na/(n_na + n_not_na)
                # If the column has greater than 20% of NA, add to drop column
                if (freq_na >= 0.2) {
                        print(var)
                        print("Less than 20%")
                        # Add the variable to drop column
                        drop_col <- c(var, drop_col) %>% unique()
                }
                if (freq_na > 0 && freq_na < 0.2) {
                        # Add variable to list of variables need na adjusting (numeric)
                        if (class(df[[var]]) %in% c('double', 'numeric', 'integer')) {
                                na_num <- c(var, na_num)
                        }
                        # Add variable to list of variables need na adjusting (categorical)
                        if (class(df[[var]]) %in% c('character', 'factor')) {
                                na_cat <- c(var, na_cat)
                        }
                        
                }
        }
        
        # Perform another iteration of dropping unusable clumns (too many NA's)
        #dependents <- names(df)[names(df) %!in% drop_col]
        
        # Correct for NA values (Categorical Variables)
        #############################
        #na_cat is null
        
        # Correct for NA values (Numerical Variables)
        #############################
        
        # Examine which variables need to be adjusted for NA values
        na_num
        summary(df %>% dplyr::select(na_num))
        str(df)
        
        # If the NA values from all variables account for a small percentage of training data, just remove those variables
        # TODO: You are removing True positives; make sure to incorporate them back into final call set
        # Determine the amount of NA across columns of interest
        n_na <- df %>% filter_at(vars(na_num), any_vars(is.na(.))) %>% nrow()
        freq_na <- n_na/(nrow(df))
        # If this is less than 2% of training data, drop the rows
        if (freq_na >= 0.02) {
                # We drop rows, they likely won't alter the model results
                df <- df %>% filter_at(vars(na_num), any_vars(!is.na(.)))
                
        }
        
        for (f in dependents) {
                # Only apply to columns that are factors
                if (is.factor(test_df[[f]])) {
                        if (sort(unique(train_df[[f]])) != sort(unique(test_df[[f]]))) {
                                print(f)
                                print('training')
                                print(unique(train_df[[f]]))
                                print('testing')
                                print(unique(test_df[[f]]))
                        }
                }
        }
        
        
        # Double check factor level length is under 53 for all factors
        lvl_cols <- c()
        # Remove columns with more than 53 categories
        for (col in dependents) {
                factor_len <- length(levels(train_df[[col]]))
                if (factor_len >= 53) {
                        lvl_cols <- append(lvl_cols, col)
                        print(factor_len)
                        print(col)
                }
        }
        
        
        ########################################################################
        
        ########################################################################
        ############## Choose the Dependent and Response Variables #############
        ########################################################################
        
        # Determine which columns should be considered in the model (response and dependents)
        response <- 'SNP_MA'
        dependents <- c('QUAL','GT','REF_AD','ALT_AD','DP','AN','AF')
        
        ########################################################################
        ######### Determine Model Algorithm to Use (Automate with Caret) #######
        ########################################################################
        
        #' ## Build Training Set
        
        ########################################################################
        ################# Split the Data #######################################
        ########################################################################
        
        # DEV:
        df <- df_mod
        
        # Perform a data split and holdout set
        
        #' ## Build Model Formulas
        
        ########################################################################
        ################# Build: Formula #######################################
        ########################################################################
        
        # Create the appropriate formula: (only important features)
        formula <- as.formula(paste(response, paste(dependents, collapse=" + "), sep=" ~ "))
        
        #' Model Formula
        print(formula)
        
        #' ## Random Forest Model Development
        #' 
        ########################################################################
        ############# Random Forest Model Development ##########################
        ########################################################################
        
        df <- df %>%
                mutate(SNP_MA = ifelse(SNP_MA == '1', 'x', 'o'))
        df$SNP_MA <- as.factor(df$SNP_MA)
        summary(df %>% dplyr::select(dependents))
        # Generate the train control
        ctrl <- trainControl(method = "repeatedcv",
                             number = 10, repeats = 1,
                             selectionFunction = "best",
                             savePredictions = TRUE,
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary)
        
        # Set up the random forest tuning grid
        grid_rf <- expand.grid(mtry = c(2))
        
        # Train the model
        set.seed(123)
        m_rf <- train(formula, 
                      data = sample_n(df, 100000),
                      method = "rf",
                      metric = "ROC",
                      trControl = ctrl,
                      tuneGrid = grid_rf)
        
        # Examine the random forest results
        m_rf
        
        # Visualize the ROC curve
        roc_rf <- roc(m_rf$pred$obs, m_rf$pred$x)
        plot(roc_rf, col = 'red', legacy.axes = TRUE)
        
        acc <- 0.9999997
        err <- 1-acc
        err*1000000000
        
        
        set.seed(123)
        summary(train_df)
        rf_model <- randomForest(formula = formula,
                                 data = train_df,
                                 na.action = na.omit,
                                 type = classification,
                                 ntree = 1000,
                                 mtry = 5)
        
        # Examine general stats about the model
        print(rf_model)
        
        # Examine Variable Importance
        #class(rf_model)
        priority_cols <- rev(importance(rf_model)[order(importance(rf_model)),])
        # Remove predictors with no affect and update formula
        priority_cols <- priority_cols[priority_cols > 0.01]
        red_deps <- names(priority_cols) 
        formula <- as.formula(paste(response, paste(red_deps, collapse=" + "), sep=" ~ "))
        length(priority_cols)
        names(priority_cols)
        varImpPlot(rf_model)
        
        # Random forest provides a built-in validation set without any extra work
        # Grab the OOB error matrix and examine
        err <- rf_model$err.rate
        head(err)
        
        # Look at final OOB error rate (last row in err matrix)
        oob_err <- err[nrow(err), "OOB"]
        print(oob_err)
        
        # Plot the model trained in the previous exercise
        plot(rf_model)
        
        # Add a legend since it doesn't have one by default
        legend(x = "topright", 
               legend = colnames(err),
               fill = 1:ncol(err))
        
        # Construct object needed to tune the model
        train_predictor_df <- train_df[ ,dependents]
        train_response_vector <- as.factor(train_df$VALIDATION)
        
        # Tune the model
        # Train the mtry parameter based on the OOB error
        trees <- c(200, 400, 1000, 2000)
        for (ntree in trees) {
                set.seed(1)
                res <- tuneRF(x = train_predictor_df,
                              y = train_response_vector,
                              ntreeTry = ntree,
                              stepFactor = .5,
                              na.action = na.omit)
        }
        
        # Examine the results
        print(res)
        
        # Manual Search (best seems to be 14 and 400 trees or 1000 trees)
        control <- trainControl(method="repeatedcv", number=10, repeats=5, search="grid")
        mtrys = c(3,4,5,6,7)
        trees <- c(1000,2000)
        modellist <- list()
        for (m in seq_along(mtrys)) {
                mtry = mtrys[m]
                print(mtry)
                tunegrid <- expand.grid(.mtry=mtry)
                for (t in seq_along(trees)) {
                        ntree = trees[t]
                        print(ntree)
                        set.seed(1)
                        fit <- caret::train(formula, data=train_df, method="rf", metric="Accuracy", 
                                            tuneGrid=tunegrid, trControl=control, ntree=ntree,
                                            na.action = na.omit)
                        key <- paste0('T_',toString(ntree),'MTRY_',toString(mtry))
                        modellist[[key]] <- fit
                } 
        }
        
        # compare results
        results <- resamples(modellist)
        summary(results)
        dotplot(results)
        
        # Landis and Koch considers 0-0.20 as slight, 0.21-0.40 as fair, 0.41-0.60 as moderate, 
        # 0.61-0.80 as substantial, and 0.81-1 as almost perfect. Fleiss considers kappas > 0.75 as
        # excellent, 0.40-0.75 as fair to good, and < 0.40 as poor.
        
        #' ##Final Generation of Model
        
        ################################################################################
        ########## Final Generation of Model ###########################################
        ################################################################################
        
        set.seed(123)
        # Generate the model with the optimal parameters
        rf_model <- randomForest(formula = formula,
                                 data = train_df,
                                 na.action = na.omit,
                                 type = classification,
                                 ntree = 1000,
                                 mtry = 5,
                                 keep.forest=TRUE,
                                 importance=TRUE)
        
        # Examine general stats about the model
        print(rf_model)
        # Random forest provides a built-in validation set without any extra work
        # Grab the OOB error matrix and examine
        err <- rf_model$err.rate
        head(err)
        # Look at final OOB error rate (last row in err matrix)
        oob_err <- err[nrow(err), "OOB"]
        print(oob_err)
        # Plot the model trained in the previous exercise
        plot(rf_model)
        
        # Add a legend since it doesn't have one by default
        legend(x = "right", 
               legend = colnames(err),
               fill = 1:ncol(err))
        
        #' ##Generate a ROC CURVE and AUC
        
        ################################################################################
        ########## Generate a ROC CURVE and AUC ########################################
        ################################################################################
        
        # http://scg.sdsu.edu/rf_r/
        # https://blog.revolutionanalytics.com/2016/11/calculating-auc.html
        # https://statquest.org/2018/12/17/roc-and-auc-in-r/ (This site is used for this script)
        
        set.seed(123)
        # Generate the model with the optimal parameters
        rf_model <- randomForest(formula = formula,
                                 data = train_df,
                                 na.action = na.omit,
                                 type = classification,
                                 ntree = 1000,
                                 mtry = 5,
                                 keep.forest=TRUE,
                                 importance=TRUE)
        
        ## pty sets the aspect ratio of the plot region. Two options:
        ##                "s" - creates a square plotting region
        ##                "m" - (the default) creates a maximal plotting region
        par(pty = "s") 
        
        # Save the plot
        #pdf('./figures/ROC_RF_all_indels_final.pdf', width = 5, height = 5)
        roc(train_df$VALIDATION, rf_model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
        #dev.off()
        
        ################################################################################
        
        #' Predict Somatic Indels from the Prediction Dataset
        
        ################################################################
        ######## Prepare the dataset for prediction ####################
        ################################################################
        
        # Extract the values in the data that have not undergone validation
        test_df2 <- dplyr::filter(df, NUM_TOOLS >= 3, VALIDATION %!in% c('TP', 'FP'))
        dim(test_df2)
        
        # Make sure to set aside the identifier columns
        restore <- c("CHROM","POS","REF","ALT",'SAMPLE')
        results_keys <- test_df[ , (names(test_df) %in% restore)]
        
        # Convert character vectors to factors
        #cvecs <- as.vector(names(test_df[, sapply(test_df, class) == 'character']))
        # Iterate through each column and transform
        #for (char_vec in cvecs) {
        #        print(char_vec)
        #print(table(factor(df[[char_vec]])))
        #        test_df[[char_vec]] <- as.factor(test_df[[char_vec]])
        #}
        #priority_cols
        
        # Use the remainder of the dataset and the model to predict
        set.seed(123)
        pred <- predict(object = rf_model,
                        newdata = test_df,
                        na.action = na.omit)
        # Apply the predicted values to the dataset
        test_df$RF_PREDICTED <- pred
        
        #' No Values were ignored (Resulting Predictions)
        # Detemrine how many values were ignored
        table(test_df$RF_PREDICTED)
        table(is.na(test_df$RF_PREDICTED))
        
        #Extract only the predicted true positives
        test_tp <- filter(test_df, RF_PREDICTED == 'TP')
        test_fp <- filter(test_df, RF_PREDICTED == 'FP')
        dim(test_tp)
        test_tp['VALIDATION'] <- test_tp['RF_PREDICTED']
        
        # Join the predicted true calls back into the original raw df
        test_tp <- as_tibble(test_tp)
        dim(test_tp)
        test_join <- test_tp %>% dplyr::select(CHROM,POS,REF,ALT,SAMPLE)
        dim(test_join)
        dim(df)
        test_orig <- inner_join(df, test_join, by= c("CHROM","POS","REF","ALT","SAMPLE"))
        dim(test_orig)
        
        # Test set
        spc <- data.frame()
        for (s in unique(test_tp$SAMPLE)){
                print(s)
                pow <- as.character(unique(filter(test_tp, SAMPLE == s) %>% dplyr::select(POWER) %>% unlist()))
                c <- nrow(filter(test_tp, SAMPLE == s))
                spc2 <- data.frame(s,pow,c)
                spc <- rbind(spc,spc2)
        }
        
        names(spc) <- c('SAMPLE','POWER','FREQ') 
        spc
        # Examine the mutation freq relationship
        ggplot(spc, aes(x=POWER, y=FREQ)) +
                geom_boxplot()
        
        # Comnine the validated and 
        test_tp <- test_tp[,names(test_tp) %!in% 'RF_PREDICTED']
        train_tp <- filter(train_df, VALIDATION == 'TP')
        
        # Combine the dataframes
        results_tp <- rbind(train_tp,test_tp)
        dim(results_tp)
        
        # No duplicated results
        #head(results_tp[duplicated(results_tp[,c('CHROM','POS','REF','ALT','SAMPLE')]),])
        
        # Extract significantly muated indel genes
        results_tp$SYMBOL <- as.character(results_tp$SYMBOL)
        sort(table(results_tp$SYMBOL))
        indel_genes <- names(table(results_tp$SYMBOL)[table(results_tp$SYMBOL) >= 1])
        # Save indel genes to file
        write.table(indel_genes, './data/indel_1_or_more_genes.txt', sep ='\t',
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        # Check for empty gene symbols
        #filter(results_tp, SYMBOL == '')
        
        # Save the results to a tab seperated file
        out_file <- "./data/indel_rf_tp_pred_final.txt"
        write.table(results_tp, out_file, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        
        # Save the results to a csv file
        out_file <- "./data/snv_indel_rf_tp_pred_final.csv"
        write.table(results_tp, out_file, sep = ",", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        
        # Important details about model
        # 52 Features used (43 features were considered significant influencers of the model)
        # Most important features were the agreement between callers, samples, altNM_T, zMQ_T, altMQ_T. refNM_N, refNM_T, NUM_TOOLS, VAF_T
        # Type of random forest: classification
        # No. of variables tried at each split: 5
        # Number of trees: 1000
        # OOB estimate of  error rate: 7.94%
        # AUC: 95.3%
        # training n: 63
        # Testing n: 4,966
        # Predicted Indels: 1,632
        # File of predcited and validated indels: "./data/indel_rf_tp_pred_final.txt"
        
        #' #### End
        
        ########################################################################
        ######################## Assess Accuracy of Model ######################
        ########################################################################
}




















#' ## Add Annotations to filtered germline SNP calls
#' ##### Annotations to Add:
#' * Number of F1 Samples with atleast 1 alternative allele called
#' * Number of F1 Samples with no alternative alleles called
#' * Number of F1 Samples successfully measured by GATK
#' 
#' * 6x7 F1 SNP calls (VSQR and hard filtered)
#'
#+ Add Annotations to filtered germline SNP calls

################################################################################
#####     Add Annotations to filtered germline SNP calls      ##################
################################################################################

# Line 6, Line 7, and 6x7 F1 germline variant (VSQR filtered) calls from WGS
#########################

# This file is about 6.5 GB and cannot be loaded into a normal R session. Instead, we will load files by chromosome.

# Extract chromosomes needed for this analysis (just those that underwent genotyping)
chromosomes <- filter(K600_df, CHROM != 'LGE22C19W28_E50C23') %>% 
        dplyr::select(CHROM) %>% unique() %>% unlist() %>% as.vector()

# 600K SNP calls previously used to train the VSQR model
# Select columns of interest from K600
K600_F1 <- K600_df %>% filter(F1_GENO_FREQ == 1) %>% dplyr::select(CHROM,POS,F1_1,F1_2,F1_3,F1_4,F1_5,F1_6,F1_GENO_FREQ)

# Examine conservation of genotypes in training set
#semi_join(K600_df, train_df, by = c('CHROM','POS')) %>% 
#        dplyr::select(F1_GENO_FREQ) %>% summary()

# Perform semi join
train_set <- semi_join(train_df, K600_F1, by = c('CHROM','POS'))
train_set <- semi_join(K600_F1, train_df,  by = c('CHROM','POS'))
View(train_set)
# Determine what percentage of the 542871 SNPs conserved across F1s were used to train
n_F1 <- dim(K600_F1)[1]
n_train <- dim(train_set)[1]
n_train/n_F1 # Almost half (48.19%)

# Generate a Test set

# Generate a test set of the 52% not in the training set
test_set <- anti_join(K600_F1, train_set, by = c('CHROM','POS'))
ncol(test_set)

# Remove any positions that agree with reference or have no SNP ID
test_tp <- test_set %>% dplyr::select('F1_1','F1_2','F1_3','F1_4','F1_5','F1_6','GG5_REF','AFFY_SNP_ID') %>% filter(!is.na(AFFY_SNP_ID))

dim(test_set)

# Adjust genotype calls to characters
for( column in names(test_tp) ){
        test_tp[[column]] <- as.character(test_tp[[column]]) 
}
dim(test_tp)
# Convert rows into lists (20 seconds)
row_list <- test_tp %>% purrr::pmap(list)
# Count the genotypes (4 minutes)
test_tp$REF_FREQ <- lapply(row_list,count_genos_f1) %>% unlist()
dim(test_tp)
test_tp <- test_tp %>% filter(REF_FREQ <= 1) %>% dplyr::select(AFFY_SNP_ID, REF_FREQ)

# Generate a revised test set
test_set <- semi_join(test_set,test_tp, by = "AFFY_SNP_ID")

# Make sure the main vcf file is bgzipped and tabix indexed

# End of pipeline
#system('bgzip ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf')
#system('tabix -p vcf ./data/vsqr/normals_recal_famprioirs_gfiltered_denovo_germline_99.9.g.vcf.gz')

# Post vsqr, post haplotyping
#system('bgzip ./data/vsqr/S1-22_germline_vsqr_99.9.g.vcf')
#system('tabix -p vcf ./data/vsqr/S1-22_germline_vsqr_99.9.g.vcf.gz')

# Post vsqr
#system('bgzip ./data/vsqr/S1-22_germline_vsqr_99.9.g.vcf')
#system('tabix -p vcf ./data/vsqr/S1-22_germline_vsqr_99.9.g.vcf.gz')

# Post Hard filters
#system('bgzip ./data/hard_filtered_variants/germline_hard_filtered_snps_lenient.g.vcf')
#system('tabix -p vcf ./data/hard_filtered_variants/germline_hard_filtered_snps_lenient.g.vcf.gz')

# Create empty chromosome to collect stats
df_out <- data.frame()
df_out
for (chrom in chromosomes) {
        print(chrom)
        print('Generating the VCF file')
        # First separate the vcf files
        vsqr_file <- paste0('./data/hard_filtered_variants/germline_hard_filtered_snps_lenient_chr',chrom,'_nohead.g.vcf')
        system_cmd <- paste0('tabix ./data/hard_filtered_variants/germline_hard_filtered_snps_lenient.g.vcf.gz ',chrom,' | grep -v "#" > ',vsqr_file)
        system(system_cmd)
        print('VCF for chromosome generated')
        
        print('Loading VCF')
        vsqr_df <- read.table.ffdf(file = vsqr_file, sep = '\t',
                                   colClasses=c("factor", "integer","factor","factor","factor","double",rep("factor",25)), VERBOSE = TRUE, first.rows = 20000 ) %>% as_tibble()
        print('VCF data Loaded')
        
        # Set the names
        names(vsqr_df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","738-0_S33","741-0_S34","756-0_S35","766-0_S36","777-0_S37","787-0_S38","788-0_S39","794-0_S40","798-0_S41","833-0_S42","834-0_S43","835-0_S44","841-0_S45","842-0_2_S28","855-0_S46","863-0_S47","884-0_S48","901-0_2_S29","906-0_S49","911-0_S30","918-0_S50","927-0_S51")
        # Remove Lines 6 and 7
        #vsqr_f1 <- vsqr_df %>% dplyr::select(-'002683_Line-6', -'002684_Line-7')
        vsqr_f1 <- vsqr_df
        # Set names
        vsqr_f1 <- vsqr_f1 %>% setNames(paste0('F1_', names(vsqr_f1)))
        names(vsqr_f1)[1:9] <- str_replace(names(vsqr_f1)[1:9], 'F1_', '')
        
        # Iterate through the columns and only select the genotype
        print('Iterating through columns to select genotypes')
        for (c in names(vsqr_f1)[10:ncol(vsqr_f1)]) {
                print(c)
                vsqr_f1 <- vsqr_f1 %>% separate(c, c(c, NA), sep = ':')
        }
        
        print("Adding Genotype Stats")
        # Transpose data frames and add annotation
        vsqr_tp <- vsqr_f1[,10:ncol(vsqr_f1)] %>% transpose()
        unlist(vsqr_f1[,10:ncol(vsqr_f1)]) %>% as.vector() %>% table()
        # Genotypes to count: ./., 0/0, 0/1, 0/2, 1/1, 1/2, 2/1, 2/2
        # Count the genotypes
        vsqr_f1$F1_n_NA <- str_count(vsqr_tp, pattern = "\\./\\.")
        vsqr_f1$F1_n_meas <- 22 - str_count(vsqr_tp, pattern = "\\./\\.")
        vsqr_f1$F1_n_00 <- str_count(vsqr_tp, pattern = '0/0')
        vsqr_f1$F1_n_01 <- str_count(vsqr_tp, pattern = '0/1')
        vsqr_f1$F1_n_02 <- str_count(vsqr_tp, pattern = '0/2')
        vsqr_f1$F1_n_11 <- str_count(vsqr_tp, pattern = '1/1')
        vsqr_f1$F1_n_12 <- str_count(vsqr_tp, pattern = '1/2')
        vsqr_f1$F1_n_21 <- str_count(vsqr_tp, pattern = '2/1')
        vsqr_f1$F1_n_22 <- str_count(vsqr_tp, pattern = '2/2')
        # Number of variants called
        # TODO: Determine cutoff for number of variants called
        vsqr_f1m <- vsqr_f1 %>% filter(F1_n_meas >= 4)
        print("Genotype Stats Finished")
        
        # Grab the chomosome from the test set
        test_chr <- test_set[test_set$CHROM == chrom,]
        train_chr <- train_set[train_set$CHROM == chrom,]
        test_in <- semi_join(test_chr, vsqr_f1m, by = c('CHROM','POS'))
        train_in <- semi_join(train_chr, vsqr_f1m, by = c('CHROM','POS'))
        
        # Collect stats for dataframe
        row_out <- list(chrom,dim(vsqr_f1m)[1],dim(train_chr)[1],dim(train_in)[1],dim(test_chr)[1],dim(test_in)[1])
        row_out <- data.frame(row_out)
        names(row_out) <- c('CHROM','VSQR_N','TRAIN_N','TRAIN_IN','TEST_N','TEST_IN')
        df_out <- rbind(df_out,row_out)
        print('Finished with Chromosome')
        print(df_out)
}



# Plot the CV 5-fold
# Also using forward stepwise selection
# Model Selection by Cross-Validation
# k-fold cross-validation.

set.seed(123)
#val <- train_df[sample(nrow(train_df), 105), ]
val<- train_set
# Collects row numbers of test set
folds = sample(rep(1:5, length = nrow(val)))
table(folds)

# Create empty lists to store results
k_list <- list()
sens_list <- list()
spec_list <- list()
auc_list <- list()

k=1
# Collect stats on k folds cross validation
for (k in 1:5) {
        set.seed(123)
        ko_data <- val[folds != k, ]
        k_data <- val[folds == k, ]
        dim(k_data)
        
        # Load the predicted data
        # Consider performing a semi join from k_data
        
        # Generate an identifier for all SNPs
        
        # Collect the FPR and TPR
        rf.pred <- prediction(pred, k_data$VALIDATION)
        # Predict the model's performance
        rf.perf <- performance(rf.pred, "tpr", "fpr")
        # Adjust the values
        rf.perf@x.values <- list(1-(rf.perf@x.values %>% unlist()))
        rf.perf@x.name <- "Specificity"
        rf.perf@y.name <- "Sensitivity"
        # Collect sensitivity and specificity
        sens <- rf.perf@y.values
        spec <- rf.perf@x.values
        # Get the AUC, TPR, and FPR
        roc_full_resolution <- roc(k_data$VALIDATION,pred)
        AUC <- as.numeric(auc(roc_full_resolution))
        
        # Store values in dataframe
        # Trail
        k_list <- append(k_list, k)
        sens_list <- append(sens_list, sens)
        spec_list <- append(spec_list, spec)
        auc_list <- append(auc_list, AUC)
}


# Create an S4 object to store data
# Set a class with setClass
setClass("ROC_Store", representation(K="list",Sensitivity="list", Specificity="list",AUC="list"),
         prototype(K=list(NA_real_),Sensitivity = list(NA_real_), Specificity = list(NA_real_), AUC=list(NA_real_)))
# Create an instance of a class with new
cv.results <- new("ROC_Store",
                  K=k_list,
                  Sensitivity = sens_list,
                  Specificity = spec_list,
                  AUC = auc_list)

k1_df <- data.frame('k'=1,'spec'=cv.results@Specificity[[1]],'sens'=cv.results@Sensitivity[[1]],auc=cv.results@AUC[[1]])
k2_df <- data.frame('k'=2,'spec'=cv.results@Specificity[[2]],'sens'=cv.results@Sensitivity[[2]],auc=cv.results@AUC[[2]])
k3_df <- data.frame('k'=3,'spec'=cv.results@Specificity[[3]],'sens'=cv.results@Sensitivity[[3]],auc=cv.results@AUC[[3]])
k4_df <- data.frame('k'=4,'spec'=cv.results@Specificity[[4]],'sens'=cv.results@Sensitivity[[4]],auc=cv.results@AUC[[4]])
k5_df <- data.frame('k'=5,'spec'=cv.results@Specificity[[5]],'sens'=cv.results@Sensitivity[[5]],auc=cv.results@AUC[[5]])
k_df <- rbind(k1_df,k2_df)
k_df$k <- as.factor(k_df$k)
# Average AUC
mean_auc <- mean(cv.results@AUC %>% unlist())
# Generate a confidence interval for the AUC
SD <- sd(cv.results@AUC %>% unlist())
SE <- SD/sqrt(length(cv.results@AUC %>% unlist()))
SE2 <- SE*2

# 0.95 +/- 0.038 95%C.I. (too high)
mean_auc

Specificity <- k3_df$spec
Sensitivity <- k3_df$sens
# Plot the ROC Curve
ggplot(k3_df,aes(x=Specificity ,y=Sensitivity, color = 'red')) +
        geom_line() +
        xlim(1,0) +
        abline()


# Save the plot
#plot(rf.perf,main="Validation set (n=106)\n(Somatic INDELs)",col=2,lwd=2,xlim=c(1,0))
