library(jsonlite)
library(ggplot2)
library(UpSetR)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(visNetwork)
library(clusterProfiler)
library(enrichplot)

# Get the present working directory
pwd <- getwd()

# Load configuration
config <- fromJSON("config.json")

# Extract paths from configuration and convert to absolute paths
outdir <- file.path(config$output_files$outdir)
graph <- file.path(outdir, config$output_files$graph)
enrichment_KEGG_results_csv <- file.path(outdir, config$output_files$enrichment_KEGG_results_csv)

setwd(graph)

# Read data from CSV file
csv_file_path <- enrichment_KEGG_results_csv
enrichment_data <- read.csv(csv_file_path)

# Get p-value threshold from command line arguments
args <- commandArgs(trailingOnly = TRUE)
pvalue_threshold <- ifelse(length(args) > 0, as.numeric(args[1]), 0.05)

# Filter data based on the provided p-value threshold
padjust_filtered_data <- subset(enrichment_data, p.adjust <= pvalue_threshold)



psdjust_enrich_res <- new("enrichResult", result = padjust_filtered_data)


psdjust_enrichres2 <- pairwise_termsim(psdjust_enrich_res) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index

png("padjust_Treeplot.png", width = 2020, height = 1080, units = "px", res = 100)

treeplot(psdjust_enrichres2, showCategory = 30, color = "p.adjust", label_format = 30)

dev.off()


