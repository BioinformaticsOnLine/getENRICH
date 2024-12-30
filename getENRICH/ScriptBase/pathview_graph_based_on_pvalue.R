library(jsonlite)
library(dplyr)
library(pathview)

# Get the present working directory
pwd <- getwd()

# Load configuration
config <- fromJSON("config.json")

# Extract paths from configuration and convert to absolute paths
outdir <- file.path(config$output_files$outdir)
pathway <- file.path(outdir, config$output_files$pathway)
enrichment_KEGG_results_csv <- file.path(outdir, config$output_files$enrichment_KEGG_results_csv)

setwd(file.path(pathway, "p_value_based_pathway"))

# Read data from CSV file
data <- read.csv(enrichment_KEGG_results_csv)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
pvalue_threshold <- ifelse(length(args) > 0, as.numeric(args[1]), 0.05)
organism <- ifelse(length(args) > 1, args[2], "ko")  # Default to "ko" if not provided

# Subset the dataframe based on the condition
df_pvalue <- data[data$pvalue < pvalue_threshold, c('ID', 'geneID')]

# Replace "ko" with the specified organism
df_pvalue$ID <- sub("^ko", organism, df_pvalue$ID)

# Loop through each row of the dataframe
for (i in 1:nrow(df_pvalue)) {
  # Extract pathway_id and ko_data for the current row
  pathway_id <- df_pvalue$ID[i]
  ko_data <- unlist(strsplit(df_pvalue$geneID[i], "/"))

  # Use tryCatch to handle errors and continue with the next iteration
  tryCatch({
    pathview(gene.data = ko_data, pathway.id = pathway_id, species = organism, plot.col.key=FALSE)
  }, error = function(e) {
    message(sprintf("Error processing pathway ID %s: %s", pathway_id, e$message))
    # Continue to the next iteration
  })
}
