library(jsonlite)
library(ggplot2)
library(UpSetR)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(visNetwork)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize) 


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



# Heatmap plot
gene_list <- strsplit(pvalue_filtered_data$geneID, "/")
gene_term_df <- data.frame(
  term = rep(pvalue_filtered_data$Description, sapply(gene_list, length)),
  gene = unlist(gene_list)
)

# Create binary presence/absence matrix for the heatmap
heatmap_data <- gene_term_df %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = gene, values_from = value, values_fill = 0)

# Debug: Inspect the heatmap data
#print("Heatmap data:")
#print(head(heatmap_data))

# Extract the matrix for heatmap plotting
heatmap_matrix <- as.matrix(heatmap_data[, -1])  # Remove `term` column
rownames(heatmap_matrix) <- heatmap_data$term

# Define color options
color_options <- list(
  BlueWhite = colorRamp2(c(0, 1), c("white", "steelblue")),
  RedYellow = colorRamp2(c(0, 1), c("white", "red")),
  GreenPurple = colorRamp2(c(0, 1), c("white", "purple"))
)

# Choose a color scheme
selected_color <- "GreenPurple"
col_fun <- color_options[[selected_color]]

# Set plot size dynamically based on data
num_pathways <- nrow(heatmap_matrix)
num_genes <- ncol(heatmap_matrix)

# Adjust heatmap width to accommodate pathway names and annotations
max_pathway_length <- max(nchar(rownames(heatmap_matrix)))  # Max length of pathway names
heatmap_height <- max(7, num_pathways * 0.4)
heatmap_width <- max(15, num_genes * 0.3, max_pathway_length * 0.2)

# Function to draw the heatmap
draw_heatmap <- function() {
  draw(
    Heatmap(
      heatmap_matrix,
      name = "Gene Presence",
      col = col_fun,                # Custom color scale
      cluster_rows = FALSE,         # Disable row clustering dendrogram
      cluster_columns = TRUE,       # Enable clustering for columns
      column_title = "Genes",       # Label for columns
      row_names_side = "left",      # Place row names on the left
      column_names_side = "top",    # Place column names on the top
      row_names_gp = gpar(fontsize = 12),  # Adjust row font size
      column_names_gp = gpar(fontsize = 9), # Adjust column font size
      top_annotation = HeatmapAnnotation(
        Pathway_Counts = anno_barplot(colSums(heatmap_matrix), gp = gpar(fill = "red"))
      ),  # Add a barplot showing the gene frequency across pathways
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 0.5, fill = fill))  # Add border and fill
      },
      border = TRUE,                  # Additional border around the heatmap
      show_heatmap_legend = FALSE     # Remove the color legend
    ),
    padding = unit(c(2, 3, 2, 3), "cm")  # Top, right, bottom, left padding
  )
}

# Save heatmap as PDF
pdf("pvalue_heatmap_plot.pdf", height = heatmap_height, width = heatmap_width,
    paper = "special", pointsize = 10)
draw_heatmap()
dev.off()

# Save heatmap as PNG
png("pvalue_heatmap_plot.png", height = heatmap_height, width = heatmap_width, units = "in", res = 300)
draw_heatmap()
dev.off()
