library(jsonlite)
library(cowplot)
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

#csv_file_path <- "enrichment_KEGG_results.csv"  # Replace with your actual file path
enrichment_data <- read.csv(csv_file_path)

# Get p-value threshold from command line arguments
args <- commandArgs(trailingOnly = TRUE)
pvalue_threshold <- ifelse(length(args) > 0, as.numeric(args[1]), 0.05)

# Filter data to include only terms with PValue <= pvalue_threshold
padjust_filtered_data <- subset(enrichment_data, p.adjust <= pvalue_threshold)

png("padjust_pubMed_trend_plot.png", width = 1920, height = 1080, units = "px", res = 100)

terms <- padjust_filtered_data$Description
p <- pmcplot(terms, 2010:2020)
p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
plot_grid(p, p2, ncol=2)

dev.off()
