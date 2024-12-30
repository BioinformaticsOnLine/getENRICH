library(jsonlite)

# Get the present working directory
pwd <- getwd()

# Load configuration
config <- fromJSON("config.json")

# Extract paths from configuration and convert to absolute paths
outdir <- file.path(config$output_files$outdir)
graph <- file.path(outdir, config$output_files$graph)
kegg_annotationTOgenes_sb_3 <- file.path(config$input_files$kegg_annotationTOgenes_sb_3)
background_genes_sb <- file.path(config$input_files$background_genes_sb_1)
genes_of_interest_sb <- file.path(config$input_files$genes_of_interest_sb_2)
enrichment_KEGG_results_csv <- file.path(outdir, config$output_files$enrichment_KEGG_results_csv)

setwd(graph)

library(dplyr)
library(tidyverse)
library(clusterProfiler)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
pvalue_threshold <- as.numeric(args[1])
organism <- args[2]

# Check if the organism is the default "ko" or a specified one
if (organism == "ko") {
  # Original process
  eggNOG_kegg <- read_tsv(kegg_annotationTOgenes_sb_3)

  background_genes <- read_tsv(background_genes_sb) %>%
    unlist() %>%
    as.vector()

  # read the gene list of interest
  interesting_set <- read_tsv(genes_of_interest_sb) %>%
    unlist() %>%
    as.vector()

  # Create a clean list of KEGG orthologs for the background
background_kegg <- eggNOG_kegg %>%
  dplyr::filter(gene %in% background_genes) %>%
  dplyr::select(term) %>%  # Select only the KEGG ortholog column
  unlist() %>%
  as.vector()

  # Create a clean list of KEGG orthologs for the genes of interest
interesting_set_kegg <- eggNOG_kegg %>%
  dplyr::filter(gene %in% interesting_set) %>%
  dplyr::select(term) %>%  # Select only the KEGG ortholog column
  unlist() %>%
  as.vector()

  enrichment_kegg <- enrichKEGG(interesting_set_kegg,
                                organism = organism,
                                keyType = "kegg",
                                pvalueCutoff = pvalue_threshold,
                                pAdjustMethod = "BH",
                                universe = background_kegg,
                                minGSSize = 10,
                                maxGSSize = 500,
                                qvalueCutoff = 0.05,
                                use_internal_data = FALSE)
} else {
  # New process for specified organism
  background <- read_tsv(background_genes_sb) %>%
    unlist() %>%
    as.vector()

  # read the gene list of interest
  foreground <- read_tsv(genes_of_interest_sb) %>%
    unlist() %>%
    as.vector()

  enrichment_kegg <- enrichKEGG(
    gene = foreground,            # Foreground Entrez Gene IDs
    organism = organism,          # Specified organism
    keyType = "ncbi-geneid",      # Entrez Gene IDs
    pvalueCutoff = pvalue_threshold,
    pAdjustMethod = "BH",
    universe = background,        # Background Entrez Gene IDs
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.05
  )
}

# Save the enrichment result
write.csv(file = paste0(enrichment_KEGG_results_csv),
          x = enrichment_kegg@result)
