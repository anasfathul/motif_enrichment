#!/usr/bin/env Rscript
library(readr)
library(tidyverse)
library(reticulate)
library(stats)
library(ggplot2)
library(pROC)

np <- import("numpy")

# Input parameteres:
# 1.) DHS_nmf file
# 2.) Metadata file
# 3.) Motif indicator file
# 4.) Output "file_name.csv" *include .csv
# 5.) Optional PRROC/PRAUC graph

args = commandArgs(trailingOnly=TRUE)

# test whether there are enough arguments provided: if not, return error message
if (length(args) < 3) {
        stop("At least something input should be provided", call.=FALSE)
}

print("Reading input matrix")
DHS_feature_nmf <- np$load(args[1])
DHS_feature_nmf_df <- data.frame(t(DHS_feature_nmf))

print("Reading Metadata File")
metadata <- read_delim(args[2], delim='\t', col_names=T)

print("Reading Motif Indicator")
motif_indicator <- read.table(args[3])
motif_indicator <- unlist(motif_indicator)

print("Running logistic regression model")
model = glm(motif_indicator ~., data=DHS_feature_nmf_df, family=binomial(link="logit"))


# Save in .csv format X1, X2, ..., X16, PRROC score, PRAUC score
write.csv(summary(model)['coefficients'],file=args[4])