#'---
#' title: "Determine Knees and Elbows for Hard Filtering of Germline SNPs and Indels"
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
#BiocManager::install("SamSPECTRAL")
#install.packages("BSgenome.Ggallus.UCSC.galGal5")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggplot2","SamSPECTRAL")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

############################################################
##### Functions ############################################
############################################################

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "steep"

## explicit gc, then execute `expr` `n` times w/o explicit gc, return timings (~"Jenny Bryan, updating work of Winston Chang")
################################################################################
benchmark <- function(n = 1, expr, envir = parent.frame()) {
        expr <- substitute(expr)
        gc()
        map(seq_len(n), ~ system.time(eval(expr, envir), gcFirst = FALSE))
}
################################################################################

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * Distributions of stats from WGS germline variant calling
#' * QD, 
#'
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# QD stats
######################

# In R, Generate the QD Plot
samples_QD = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_QD.txt")

# Place data in a vector
QD <- as.vector(samples_QD$V1)
# Calculate and plot the density
QD_den <- density(QD)
length(QD)
length(QD_den$y)
str(QD_den)

plot(QD_den)

# Function to get elbow/knee point indicices; AKA takes second derivative to find maximum curvature on the curve
# Function from Sandipan Dey: https://stackoverflow.com/questions/41518870/finding-the-elbow-knee-in-a-curve
################################################################################
get.elbow.points.indices <- function(x, y, threshold) {
        d1 <- diff(y) / diff(x) # first derivative
        d2 <- diff(d1) / diff(x[-1]) # second derivative
        indices <- which(abs(d2) > threshold)  
        indices
}
################################################################################

# first approximate the function, since we have only a few points
x <- QD_den$x
y <- QD_den$y
ap <- approx(x, y, n=1000, yleft=min(y), yright=max(y))
x <- ap$x
y <- ap$y

indices <- get.elbow.points.indices(x, y, 1e4) # threshold for huge jump = 1e4
x[indices]
#[1] 6.612851 # there is one such point
plot(x, y, pch=19)
points(x[indices], y[indices], pch=19, col='red')

# first approximate the function, since we have only a few points
ap <- approx(x, y, n=1000, yleft=min(y), yright=max(y))
x <- ap$x
y <- ap$y

indices <- get.elbow.points.indices(x, y, 1e4) # threshold for huge jump = 1e4
x[indices]
#[1] 6.612851 # there is one such point
plot(x, y, pch=19)
points(x[indices], y[indices], pch=19, col='red')


## Data
values <- rep(1,times=10)
values <- c(values,(10:0)/10)

## Looks like knee point:
plot(values)

## Find the knee point:
detected <- kneepointDetection(vect=values, PlotFlag=FALSE)
detected <- kneepointDetection(vect=QD_den$y, PlotFlag=FALSE)

print(detected$MinIndex)
plot(QD_den)

# FS
##########################
# Load in the FS data
samples_FS = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_FS.txt")

# Place data in a vector
FS <- as.vector(samples_FS$V1)
# Calculate and plot the density
FS_den <- density(FS)
plot(FS_den)
?plot
# Find the knee point
detected <- kneepointDetection(vect=FS_den$y, PlotFlag=TRUE)

?kneepointDetection

print(detected$MinIndex)
plot(FS_den) + xlim(0, 10)
qplot(V1, data = samples_FS, geom = "density") + xlim(0, 10.0)
qplot(V1, data = samples_FS, geom = "density") + xlim(0, 5.0)


# In R, Generate the QD Plot
samples_QD = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_QD.txt")
qplot(V1, data = samples_QD, geom = "density") + xlab("QD") + ylab("Density")

ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_QD_density_snps.png")

# In R, Generate the FS Plot (x-axis is logged)
samples_FS = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_FS.txt")
qplot(V1, data = samples_FS, geom = "density", log="x") + xlab("FS") + ylab("Density")

ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_FS_density_snps.png")
# Cutoff right after right peak

# In R, Generate the MQ Plot
samples_MQ = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_MQ.txt")
qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_snps.png")

# Zoom in on the MQ peak at 40
qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(39.0, 41.0)
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_40_snps.png")

# Zoom in on the MQ peak at 60
qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(59.0, 61.0)
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_60_snps.png")
# Discard anything that is not 60

# In R, Generate the MQRankSum Plot
samples_MQRankSum = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_MQRankSum.txt")
qplot(V1, data = samples_MQRankSum, geom = "density") + xlab("MQRankSum") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQRankSum_density_snps.png")
# Remove anything less than -2

# In R, Generate the ReadPosRankSum Plot
samples_ReadPosRankSum = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_ReadPosRankSum.txt")

# Place data in a vector
RPRS <- as.vector(samples_ReadPosRankSum$V1)
# Calculate and plot the density
RPRS_den <- density(RPRS)
plot(RPRS_den)



qplot(V1, data = samples_ReadPosRankSum, geom = "density") + xlab("ReadPosRankSum") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_ReadPosRankSum_density_snps.png")
