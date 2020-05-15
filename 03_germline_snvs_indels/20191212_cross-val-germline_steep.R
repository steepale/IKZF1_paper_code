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

# Function to load germline vcfs
################################################################################
load_normal_vcfs <- function(x, extract_gt = TRUE) {
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
cat_mode <- function(x) {
        uniqx <- unique(na.omit(x))
        uniqx[which.max(tabulate(match(x, uniqx)))]
}
################################################################################

# Generate a frequency of genotypes
################################################################################
count_genos_f1 <- function(x) {
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

################################################################################
#####     Cross Validation Tutorial      #######################################
################################################################################

# Bootstrap validation and cross validation results in similar results in practice

library(caret)
# Set sedd for reproducibility
set.seed(42)
mtcars
?train
# Fit different models and perform cross-validation (10 fold)
model <- train(
        mpg ~ hp, mtcars,
        method = 'lm',
        trControl = trainControl(
                method = "cv",
                number = 10,
                verboseIter = TRUE
        )
)






#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Line 6, Line 7, and 6x7 F1 Genotype calls from 600K SNP Arrays with genotype similarities
#' * 600K Genotypes: training set for the VSQR model
#' * Line 6, Line 7, and 6x7 F1 germline variant (VSQR filtered) calls from WGS
#' * Line 6, Line 7, and 6x7 F1 germline variant (hard filtered) calls from WGS
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

# Adjust the CHROM column names
names(K600_df)[2] <- "CHROM"

# Filter out any non-informative genotypes
K600_df <- K600_df %>% filter(CHROM != '---' | !is.na(POS))

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

#' ## Perform liftover from Galgal4 to Galgal5 on Arkansas SNPs
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
K600_df_gr <- K600_df_gr %>% filter(CHROM != "LGE22C19W28_E50C23")

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

#' ## Add Annotations to filtered germline SNP calls
#' ##### Annotations to Add:
#' * Number of F1 Samples with atleast 1 alternative allele called
#' * Number of F1 Samples with no alternative alleles called
#' * Number of F1 Samples successfully measured by GATK
#' 
#' * 6x7 F1 SNP calls (VSQR and hard filtered)
#'
#+ Load the Data

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
