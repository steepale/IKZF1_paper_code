#'---
#' title: "Predictive Model of Germline SNPs: Transform, Split, Choose Model"
#' author: Alec Steep
#' date: "20200203"
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
#' * Transform data as needed (center,scale,box-cox,remove zero variance predictors)
#' * Choose predictors
#' * Split data accordingly
#' * Choose the model algorithm

#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/MDV_Project/germline_snps_indels'

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("ff")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","ROCR","ff","GenomicRanges","BSgenome","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal5","rtracklayer", "caret","pROC","modelr","ggplot2","e1071","doMC","glmnet", "lattice","MASS","pamr","pls","sparseLDA","lubridate","reshape2","kernlab", "klaR","latticeExtra","earth","partykit","gtools")
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
#' * Processed data for model training
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Data previously processed in: 20191215_data-prep-germline_steep.R

# Processed data
######################
proc_file <- paste0(WD,'/data/processed-data_steep.txt')
col_class <- c("factor", "integer",rep("factor",2),"double","factor",rep("integer",4),
               rep("factor",2), rep('integer',4),rep("double",7),"integer",rep("double",6),
               "integer",rep("factor",6))
df <- read.table(file = proc_file, sep = '\t', header=TRUE) %>% as_tibble()

# Adjust the values of the response variable
df <- df %>%
        mutate(SNP_MA = ifelse(SNP_MA == '1', 'Y', 'N'))
df$SNP_MA <- as.factor(df$SNP_MA)

################################################################################
############## Choose the Dependent and Response Variables #####################
################################################################################

# Determine which columns should be considered in the model (response and dependents)
response <- 'SNP_MA'
dependents_grouped <- c('QUAL','GT','REF_AD','ALT_AD','DP','AN','AF')

################################################################################
############## Create Dummy Variables for Categorical Dependents ###############
################################################################################

# Apply dummary variables
dmy <- dummyVars(" ~ GT ", data = df)
dum_df <- as_tibble(predict(dmy, newdata = df))
# The results are returned as numeric factors, make them factors
dum_df$`GT../.` <- as.factor(dum_df$`GT../.`)
dum_df$`GT.0/0` <- as.factor(dum_df$`GT.0/0`)
dum_df$`GT.0/1` <- as.factor(dum_df$`GT.0/1`)
dum_df$`GT.1/1` <- as.factor(dum_df$`GT.1/1`)
# Bind back into main dataframe
df <- cbind(df, dum_df)

# Create a list of dependent variables (ungrouped as opposed to grouped)
dependents_ungrouped <- c('QUAL','GT../.', 'GT.0/0' , 'GT.0/1', 'GT.1/1' ,'REF_AD','ALT_AD','DP','AN','AF')

# Create the grouped and ungrouped dataframes
df_grouped <- df %>% dplyr::select(all_of(dependents_grouped), all_of(response))
df_ungrouped <- df %>% dplyr::select(all_of(dependents_ungrouped), all_of(response))

################################################################################
############## Visualize and Transform Certain Variables #######################
################################################################################

# Visualize before transformations
#####################
hist(df_grouped$REF_AD, breaks = 100)
hist(df_grouped$ALT_AD, breaks = 100)

# Grouped Predictors
#####################
## Use caret's preProcess function to transform for skewness
# Center and Scale, BoxCox, Remove zero variance predictors
pp_grouped <- preProcess(df_grouped[,dependents_grouped], 
                         method = c("center", "scale", "BoxCox", "nzv"))

## Apply the transformations
df_grouped <- predict(pp_grouped, newdata = df_grouped) %>% as_tibble()
df_grouped$SNP_MA <- as.factor(df_grouped$SNP_MA)

# Ungrouped Predictors
#####################

## Use caret's preProcess function to transform for skewness
# Center and Scale, BoxCox, Remove zero variance predictors
pp_ungrouped <- preProcess(df_ungrouped[,dependents_ungrouped], 
                           method = c("center", "scale", "BoxCox", "nzv"))
#pp_ungrouped <- preProcess(df_ungrouped[,dependents_ungrouped], 
#                           method = c("nzv"))

## Apply the transformations
df_ungrouped <- predict(pp_ungrouped, newdata = df_ungrouped) %>% as_tibble()
df_ungrouped$SNP_MA <- as.factor(df_ungrouped$SNP_MA)

# Visualize after transformations
#####################
hist(df_grouped$REF_AD, breaks = 100)
hist(df_grouped$ALT_AD, breaks = 100)

# Update the dependents
dependents_grouped <- names(df_grouped)[names(df_grouped) %!in% c("SNP_MA")]
dependents_ungrouped <- names(df_ungrouped)[names(df_ungrouped) %!in% c("SNP_MA")]

################################################################################
############## Test Predictors for High Correlation ############################
################################################################################

# Grouped
#####################
predCorr <- cor(df_grouped[,c("REF_AD", "ALT_AD", "DP")])
highCorr <- findCorrelation(predCorr, .99)
# TODO: Clean up this command
#df_grouped <- df_grouped[-highCorr]

# Ungrouped
#####################
predCorr <- cor(df_ungrouped[,c("REF_AD", "ALT_AD", "DP")])
highCorr <- findCorrelation(predCorr, .99)
# TODO: Clean up this command
#df_ungrouped <- df_ungrouped[-highCorr]

################################################################################
#####     Split the Data (Train/Test/Holdouts/Validation)      #################
################################################################################

# caret provides the createDataPartition() function that creates partitions based on stratified holdout sampling

# Grouped
########################
set.seed(123)
train_rows <- createDataPartition(df_grouped$SNP_MA, p = 0.70, list = FALSE) %>% as.vector()
train_grouped <- df_grouped[train_rows, ]
test_grouped <- df_grouped[-train_rows, ]

# Ungrouped
#########################
set.seed(123)
train_rows <- createDataPartition(df_ungrouped$SNP_MA, p = 0.70, list = FALSE) %>% as.vector()
train_ungrouped <- df_ungrouped[train_rows, ]
test_ungrouped <- df_ungrouped[-train_rows, ]

dim(train_ungrouped)
dim(test_ungrouped)



########################################################################
######### Save the Training/Set Data for Later Use #####################
########################################################################
# Create a train and test set to save
train_set <- df[train_rows, ]
test_set <- df[-train_rows, ]
# Save the train file
train_file <- paste0(WD,'/data/20200126_germline-snp-ml-train-set-all-ann_steep.txt')
write.table(train_set, file=train_file, sep = '\t', quote = FALSE, row.names = FALSE)
# Save the test file
test_file <- paste0(WD,'/data/20200126_germline-snp-ml-test-set-all-ann_steep.txt')
write.table(test_set, file=test_file, sep = '\t', quote = FALSE, row.names = FALSE)

########################################################################
######### Determine Model Algorithm to Use (Automate with Caret) #######
########################################################################

#' ## Build Model Formula(s)

# Create the appropriate formula: (only important features)
#formula <- as.formula(paste(response, paste(dependents, collapse=" + "), sep=" ~ "))

#' Model Formula
#print(formula)

#' ## Random Forest Model Development
#' 
################################################################################
############# Random Forest Model Development ##################################
################################################################################

# Grouped
####################
#train_grouped <- sample_n(train_grouped, 100000)
train_grouped$SNP_MA <- as.factor(train_grouped$SNP_MA)
train_grouped <- as.data.frame(train_grouped)

# Generate the train control
ctrl <- trainControl(method = "repeatedcv",
                     number = 10, repeats = 10,
                     selectionFunction = "best",
                     savePredictions = TRUE,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

# Set up the random forest tuning grid
grid_rf <- expand.grid(mtry = c(1))

# Train the model
set.seed(123)
start_time <- Sys.time()
rf_fit_grouped <- train(x = train_grouped[,dependents_grouped], 
                        y = train_grouped[,response],
                        method = "rf",
                        metric = "ROC",
                        trControl = ctrl,
                        tuneGrid = grid_rf)
end_time <- Sys.time()
end_time - start_time

# Examine the random forest results
rf_fit_grouped$bestTune

# Predict how well the model performed on the test data
test_grouped$pred <- predict(rf_fit_grouped, test_grouped)
test_grouped$pred_prob <- predict(rf_fit_grouped, test_grouped, type = "prob")

# Determine how model performs best on test set
test_grouped$pred <- as.factor(test_grouped$pred)
test_grouped$SNP_MA <- as.factor(test_grouped$SNP_MA)
confusionMatrix(test_grouped$pred,
                test_grouped$SNP_MA, positive = "Y")
table(actual = test_grouped$SNP_MA, predicted = test_grouped$pred)


# Ungrouped
####################

# CAUTION: SMALL SAMPLE
#train_ungrouped <- sample_n(train_ungrouped, 100000)
train_ungrouped$SNP_MA <- as.factor(train_ungrouped$SNP_MA)
train_ungrouped <- as.data.frame(train_ungrouped)

# Generate the train control
ctrl <- trainControl(method = "repeatedcv",
                     number = 10, repeats = 10,
                     selectionFunction = "best",
                     savePredictions = TRUE,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

# Set up the random forest tuning grid
grid_rf <- expand.grid(mtry = c(1))

# Train the model
set.seed(123)
start_time <- Sys.time()
rf_fit_ungrouped <- train(x = train_ungrouped[,dependents_ungrouped], 
                          y = train_ungrouped[,response],
                          method = "rf",
                          metric = "ROC",
                          trControl = ctrl,
                          tuneGrid = grid_rf)
end_time <- Sys.time()
end_time - start_time


# Examine the random forest results
rf_fit_ungrouped$bestTune

# Predict how well the model performed on the test data
test_ungrouped$pred <- predict(rf_fit_ungrouped, test_ungrouped)
test_ungrouped$pred_prob <- predict(rf_fit_ungrouped, test_ungrouped, type = "prob")

# Determine how model performs best on test set
test_ungrouped$pred <- as.factor(test_ungrouped$pred)
test_ungrouped$SNP_MA <- as.factor(test_ungrouped$SNP_MA)
confusionMatrix(test_ungrouped$pred,
                test_ungrouped$SNP_MA, positive = "Y")
table(actual = test_ungrouped$SNP_MA, predicted = test_ungrouped$pred)

# Examine Importance of predictors
rf_imp_ungrouped <- varImp(rf_fit_ungrouped, scale = FALSE)
plot(rf_imp_ungrouped)

################################################################################
###################### CART Model Development ##################################
################################################################################

# Generate the train control
ctrl <- trainControl(method = "repeatedcv",
                     number = 10, repeats = 10,
                     selectionFunction = "best",
                     savePredictions = TRUE,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

# Ungrouped
#############################
# Train the model
train_ungrouped <- as.data.frame(train_ungrouped)
set.seed(123)
start_time <- Sys.time()
rpart_fit_ungrouped <- train(x = train_ungrouped[,dependents_ungrouped], 
                             y = train_ungrouped[,response],
                             method = "rpart",
                             tuneLength = 30,
                             metric = "ROC",
                             trControl = ctrl)
end_time <- Sys.time()
end_time - start_time
train_ungrouped <- as_tibble(train_ungrouped)

# Grouped
#############################
# Train the model
set.seed(123)
start_time <- Sys.time()
rpart_fit_grouped <- train(x = train_grouped[,dependents_grouped], 
                           y = train_grouped[,response],
                           method = "rpart",
                           tuneLength = 30,
                           metric = "ROC",
                           trControl = ctrl)
end_time <- Sys.time()
end_time - start_time

# Examine the best complexity parameter
rpart_fit$bestTune

# Proceed with the ungrouped analysis: model performance is identical but easier to interpret
rpart_fit <- rpart_fit_ungrouped
# Save the model to a file
#saveRDS(rpart_fit, paste0("./models/",date,"_cart-model-germline-snps_steep.rds"))

# To load the model
rpart_fit <- readRDS(paste0(WD,"/models/20200203_cart-model-germline-snps_steep.rds"))

# Examine the best complexity parameter
rpart_fit$bestTune

################################################################################
######################## Model Evaluation ######################################
################################################################################

# The trainging statistics of the model were manually saved to: ./data/20200203_cart-model-train-stats_steep.txt

# Predict how well the model performed on the test data
test_ungrouped$pred <- predict(rpart_fit, test_ungrouped)
test_ungrouped$pred_prob <- predict(rpart_fit, test_ungrouped, type = "prob")

# Determine how model performs best on test set
test_ungrouped$pred <- as.factor(test_ungrouped$pred)
test_ungrouped$SNP_MA <- as.factor(test_ungrouped$SNP_MA)
rpart_fit_cm <- confusionMatrix(test_ungrouped$pred,
                                test_ungrouped$SNP_MA, positive = "Y")

# Save the confusion matrix to a file
cm_out <- data.frame(cbind(t(rpart_fit_cm$overall),t(rpart_fit_cm$byClass)))
#write.table(cm_out,
#            file=paste0("./data/20200203_test-confusion-matrix-cart-germline-snps_steep.txt"),
#            sep = '\t', quote = FALSE, row.names = FALSE)

# Visualize the tree (also save to file)
#pdf(paste0('./plots/',date,'_cart-model-vis-transformed_steep.pdf'),width=26,height=14)
plot(as.party(rpart_fit$finalModel))
#dev.off()

# Precision and recall
prec <- posPredValue(test_ungrouped$pred,
                     test_ungrouped$SNP_MA, positive = "Y")
rec <- sensitivity(test_ungrouped$pred,
                   test_ungrouped$SNP_MA, positive = "Y")
# The F-score (harmonic mean)
f <- (2 * prec * rec) / (prec + rec)
f

# Plot a roc curve
roc_rpart_train <- roc(rpart_fit$pred$obs, rpart_fit$pred$Y)
roc_rpart_test <- roc(test_ungrouped$SNP_MA, test_ungrouped$pred_prob$Y)

#pdf(paste0('./plots/',date,'_cart-roc-test_steep.pdf'),width=6,height=6)
## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region
##                "m" - (the default) creates a maximal plotting region
par(pty = "s") 
#plot(roc_rpart_train, col = 'blue', legacy.axes = TRUE, print.auc=TRUE)
plot(roc_rpart_test, col = 'red', legacy.axes = TRUE, print.auc = TRUE)
#auc(roc_rpart)
#dev.off()

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

#'The F-measure (F1-score, F-score):
#'A measure of model performance that combines precision and recall into a single number (harmonic mean).
#'
#'F-score = (2 x precision x recall) / (recall + precision)
#'
#'~Lantz, B. Machine Learning with R. (Packt, 2019)


################################################################################
######################## Summary of Analysis ###################################
################################################################################

#' Methodology:
#' * 11 Predictors were originally considered including numerical and categorical
#' * Categorical predictors were coerced into dummy variables
#' * Numeric predictors were centered, scaled
#' * Box-cox transformations were performed on 2 numeric predictors
#' * 3 Numeric predictors were removed for zero variance
#' * Final model was built with 7 predictors (1 predictors as 4 dummy variables + 3 numeric predictors)
#' * SNPs from 22 F1 normal tissues were concatenated. Data was fairly homogenous and was particioned using the createDataPartition function in caret. 70% (n=199,094) of the data was used to generate a training set, while 30% (n=85,325) were used as a test set.
#' * Random forest and CART tress were examined for analysis. Surprisingly, CART performance demonstrated equivalent performance to random forest. CART was therefore choosen for its simplicity, interpritability, and computational efficiency.
#' CART complexity parameter was 0, unpruned trees. Likely because of low predictor amount.
#' Models were evaluated on training set via 10x repeated 10-fold cross validation.

#' Summary Files and Plots:
#' Final CART Model:
#' * ./data/20200203_cart-model-germline-snps_steep.rds
#' Training Stats of CART Model:
#' * ./data/20200203_cart-model-train-stats_steep.txt
#' Test Set Confusion Matrix:
#' * ./data/20200203_test-confusion-matrix-cart-germline-snps_steep.txt
#' Model Visualizations (both transformed and untransformed)
#' * ./data/20200203_cart-model-vis-transformed_steep.pdf
#' * ./data/20200203_cart-model-vis-untransformed_steep.pdf
#' ROC Curve with AUC for Test Set:
#' * ./plots/20200203_cart-roc-test_steep.pdf






