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
pvalue_filtered_data <- subset(enrichment_data, pvalue <= pvalue_threshold)

# Prepare data for UpSet plot
gene_list <- strsplit(pvalue_filtered_data$geneID, "/")
gene_term_df <- data.frame(
  term = rep(pvalue_filtered_data$Description, sapply(gene_list, length)),
  gene = unlist(gene_list)
)

# Create a binary matrix for the UpSet plot
upset_data <- gene_term_df %>%
  mutate(value = 1) %>%
  spread(term, value, fill = 0)

# Determine the number of terms for scaling text size and plot dimensions
num_terms <- nrow(pvalue_filtered_data)

# Scale text sizes based on the number of terms
base_size <- 2  # Base size for text
text_scale_factor <- max(1, min(3, num_terms / 10))  # Scale factor, adjust as necessary
text_sizes <- c(base_size * text_scale_factor, base_size * text_scale_factor,
                base_size * text_scale_factor, base_size * text_scale_factor,
                base_size * text_scale_factor * 0.8, base_size * text_scale_factor)

# Calculate plot dimensions based on the number of terms
base_width <- 3000  # Base width
base_height <- 2800  # Base height
width_factor <- max(1, num_terms / 5)  # Adjust width based on number of terms
height_factor <- max(1, num_terms / 5)  # Adjust height based on number of terms
plot_width <- base_width * width_factor
plot_height <- base_height * height_factor

# Save the UpSet plot as a PNG file with custom dimensions and high resolution
png("pvalue_upset_plot.png", width = plot_width, height = plot_height, res = 300)
upset(
  upset_data,
  sets = colnames(upset_data)[-1],
  order.by = "freq",
  sets.bar.color = "skyblue",
  matrix.color = "red",
  main.bar.color = "blue",
  text.scale = text_sizes,  # Use the scaled text sizes
  keep.order = TRUE,  # Keep the order of sets as in the data
  point.size = 4,  # Adjust the size of the points in the matrix
  line.size = 0.5  # Adjust the size of the lines in the matrix
)
dev.off()

# Save the UpSet plot as a PDF file with custom dimensions
pdf("pvalue_upset_plot.pdf", width = plot_width / 300, height = plot_height / 300)  # Convert dimensions to inches
upset(
  upset_data,
  sets = colnames(upset_data)[-1],
  order.by = "freq",
  sets.bar.color = "skyblue",
  matrix.color = "red",
  main.bar.color = "blue",
  text.scale = text_sizes,  # Use the scaled text sizes
  keep.order = TRUE,  # Keep the order of sets as in the data
  point.size = 4,  # Adjust the size of the points in the matrix
  line.size = 0.5  # Adjust the size of the lines in the matrix
)
dev.off()
