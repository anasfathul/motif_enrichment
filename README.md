# Motif Enrichment R Script

This R script, motif_enrichment.R, performs motif enrichment analysis using logistic regression model. It takes input parameters and generates output files with results.

- Prerequisites
Before running the script, make sure you have the following R packages installed:
readr
tidyverse
reticulate
stats
PRROC

- Input Parameters
The script expects the following input parameters:

1. DHS_nmf file: The file containing the DHS_nmf data.
2. Motif indicator file: The file containing the motif indicator data.
3. Output prefix: The prefix of the output filenames.
4. Optional PRROC/PRAUC graph: An optional parameter for PRROC/PRAUC graph.
5. Metadata file: The file containing metadata information.

- Running the Script
To run the script, use the following command:
Rscript motif_enrichment.R arg1 arg2 arg3 output.csv

- Output
The script generates the following output files:

<prefix>.prroc.tsv: A CSV file containing the results of PRROC analysis.
<prefix>.coeff.tsv: A CSV file containing the coefficients from the logistic regression model.

- Description
The script performs the following steps:

Reads the input matrix and metadata file.
Reads the motif indicator data.
Runs a logistic regression model using the motif indicator and input matrix data.
Calculates the coefficients and saves them in a dataframe.
Calculates the area under the ROC curve (AUC) and area under the precision-recall curve (PR-AUC).
Saves the results in separate CSV files