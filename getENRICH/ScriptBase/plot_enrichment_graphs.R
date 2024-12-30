library(jsonlite)
library(ggplot2)
library(UpSetR)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(visNetwork)
library(clusterProfiler)
library(enrichplot)
library(plotly)

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
padjust_filtered_data <- subset(enrichment_data, p.adjust <= pvalue_threshold)

# Calculate the number of unique terms in the Description column in p-value data
pvalue_num_terms <- length(unique(pvalue_filtered_data$Description))

# Calculate the number of unique terms in the Description column in p.adjust data
padjust_num_terms <- length(unique(padjust_filtered_data$Description))

# Set dynamic text size based on the number of terms in p-value plots
if (pvalue_num_terms > 50) {
  pvalue_text_size <- 6
} else if (pvalue_num_terms > 30) {
  pvalue_text_size <- 10
} else if (pvalue_num_terms > 20) {
  pvalue_text_size <- 12
} else {
  pvalue_text_size <- 14
}

# Set dynamic text size based on the number of terms in p.adjust plots
if (padjust_num_terms > 50) {
  padjust_text_size <- 6
} else if (padjust_num_terms > 30) {
  padjust_text_size <- 10
} else if (padjust_num_terms > 20) {
  padjust_text_size <- 12
} else {
  padjust_text_size <- 14
}




# Save the bar plot as a PNG file with custom dimensions
png("p_value_bar_plot.png", width = 1200, height = 1900, res = 150)
ggplot(pvalue_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P-Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()






# Save the bar plot as a PDF file with custom dimensions
pdf("p_value_bar_plot.pdf", width = 8, height = 12)
ggplot(pvalue_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P-Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()




# Save the bar plot as a PNG and PDF file with custom dimensions

# Save the bar plot as a PNG file with custom dimensions
png("p_adjust_bar_plot.png", width = 1200, height = 1900, res = 150)
ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()  # Close the PNG device

# Save the bar plot as a PDF file with custom dimensions
pdf("p_adjust_bar_plot.pdf", width = 8, height = 12)
ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()  # Close the PDF device

# Interactive p-value bar plot
p <- ggplot(pvalue_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P-Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = pvalue_text_size),
    axis.text.x = element_text(size = 5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_value_bar_plot.html")

# Save the bar plot as a PNG and PDF file with custom dimensions

# Save the bar plot as a PNG file with custom dimensions

png("p_adjust_bar_plot.png", width = 1200, height = 1900, res = 150)
ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()



# Save the bar plot as a PDF file with custom dimensions

pdf("p_adjust_bar_plot.pdf", width = 8, height = 12)
ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive p.adjust bar plot
p <- ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "Gene Count",
       y = "Enriched Terms",
       fill = "p.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = padjust_text_size),
    axis.text.x = element_text(size = 5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_adjust_bar_plot.html")





# Save the bar plot as a PNG and PDF file with custom dimensions
# Save the bar plot as a PNG file with custom dimensions

png("p_value_bar_plot_qscore.png", width = 1200, height = 1900, res = 150)
ggplot(pvalue_filtered_data, aes(x = qvalue, y = reorder(Description, qvalue), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "qscore",
       y = "Enriched Terms",
       fill = "P-Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()





# Save the bar plot as a PDF file with custom dimensions

pdf("p_value_bar_plot_qscore.pdf", width = 8, height = 12)
ggplot(pvalue_filtered_data, aes(x = qvalue, y = reorder(Description, qvalue), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "qscore",
       y = "Enriched Terms",
       fill = "P-Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive p-value bar plot qscore
p <- ggplot(padjust_filtered_data, aes(x = qvalue, y = reorder(Description, qvalue), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "qvalue",
       y = "Enriched Terms",
       fill = "P-Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = pvalue_text_size),
    axis.text.x = element_text(size = 5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_value_bar_plot_qscore.html")







# Save the bar plot as a PNG and PDF file with custom dimensions
# Save the bar plot as a PNG file with custom dimensions

png("p_adjust_bar_plot_qscore.png", width = 1200, height = 1900, res = 150)
ggplot(padjust_filtered_data, aes(x = qvalue, y = reorder(Description, qvalue), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "qscore",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()




# Save the bar plot as a PDF file with custom dimensions

pdf("p_adjust_bar_plot_qscore.pdf", width = 8, height = 12)
ggplot(padjust_filtered_data, aes(x = qvalue, y = reorder(Description, qvalue), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "qscore",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive padjust bar plot qscore
p <- ggplot(padjust_filtered_data, aes(x = qvalue, y = reorder(Description, qvalue), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment Analysis Bar Plot",
       x = "qvalue",
       y = "Enriched Terms",
       fill = "P.adjust") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = padjust_text_size),
    axis.text.x = element_text(size = 5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_adjust_bar_plot_qscore.html")







# Save the dot plot as a PNG and PDF file with custom dimensions
# Save the bar plot as a PNG file with custom dimensions

png("p_value_dot_plot.png", width = 1200, height = 1900, res = 150)
ggplot(pvalue_filtered_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = pvalue)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(5, 10)) +
  labs(title = "Enrichment Analysis Dot Plot",
       x = "Gene Ratio",
       y = "Enriched Pathway",
       color = "P-Value",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()






# Save the bar plot as a PDF file with custom dimensions

pdf("p_value_dot_plot.pdf", width = 8, height = 12)
ggplot(pvalue_filtered_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = pvalue)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(5, 10)) +
  labs(title = "Enrichment Analysis Dot Plot",
       x = "Gene Ratio",
       y = "Enriched Pathway",
       color = "P-Value",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive p-value dot plot
p <- ggplot(pvalue_filtered_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = pvalue)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(5, 10)) +
  labs(title = "Enrichment Analysis Dot Plot",
       x = "Gene Ratio",
       y = "Enriched Pathway",
       color = "P-Value",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_value_dot_plot.html")







# Save the dot plot as a PNG and PDF file with custom dimensions
# Save the bar plot as a PNG file with custom dimensions

png("p_adjust_dot_plot.png", width = 1200, height = 1900, res = 150)
ggplot(padjust_filtered_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(5, 10)) +
  labs(title = "Enrichment Analysis Dot Plot",
       x = "Gene Ratio",
       y = "Enriched Pathway",
       color = "P.adjust",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()



# Save the bar plot as a PDF file with custom dimensions

pdf("p_adjust_dot_plot.pdf", width = 8, height = 12)
ggplot(padjust_filtered_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(5, 10)) +
  labs(title = "Enrichment Analysis Dot Plot",
       x = "Gene Ratio",
       y = "Enriched Pathway",
       color = "P.adjust",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive padjust dot plot
p <- ggplot(padjust_filtered_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(5, 10)) +
  labs(title = "Enrichment Analysis Dot Plot",
       x = "Gene Ratio",
       y = "Enriched Pathway",
       color = "P-adjust",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_adjust_dot_plot.html")






# Save the Lollipop plot as a PNG and PDF file with custom dimensions
# Save the bar plot as a PNG file with custom dimensions

png("padjust_lollipop_plot.png", width = 1200, height = 1900, res = 150)
ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(xend = 0, yend = Description), size = 0.7) +  # Increase line size
  geom_point(aes(size = Count, color = p.adjust)) +  # Increase dot size
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +  # Set minimum and maximum size for dots
  labs(title = "Enrichment Analysis Lollipop Plot",
       x = "Gene Count",
       y = "Enriched Pathway",
       color = "p.adjust",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()




# Save the bar plot as a PDF file with custom dimensions

pdf("padjust_lollipop_plot.pdf", width = 8, height = 12)
ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(xend = 0, yend = Description), size = 0.7) +  # Increase line size
  geom_point(aes(size = Count, color = p.adjust)) +  # Increase dot size
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +  # Set minimum and maximum size for dots
  labs(title = "Enrichment Analysis Lollipop Plot",
       x = "Gene Count",
       y = "Enriched Pathway",
       color = "p.adjust",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive padjust lollipop plot
p <- ggplot(padjust_filtered_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(xend = 0, yend = Description), size = 0.7) +  # Increase line size
  geom_point(aes(size = Count, color = p.adjust)) +  # Increase dot size
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 7)) +  # Set minimum and maximum size for dots
  labs(title = "Enrichment Analysis Lollipop Plot",
       x = "Gene Count",
       y = "Enriched Pathway",
       color = "p.adjust",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = padjust_text_size),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "padjust_lollipop_plot.html")







# Save the Lollipop plot as a PNG and PDF file with custom dimensions
# Save the bar plot as a PNG file with custom dimensions

png("p_value_lollipop_plot.png", width = 1200, height = 1900, res = 150)
ggplot(pvalue_filtered_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(xend = 0, yend = Description), size = 0.7) +  # Increase line size
  geom_point(aes(size = Count, color = pvalue)) +  # Increase dot size
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +  # Set minimum and maximum size for dots
  labs(title = "Enrichment Analysis Lollipop Plot",
       x = "Gene Count",
       y = "Enriched Pathway",
       color = "P-Value",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()






# Save the bar plot as a PDF file with custom dimensions

pdf("p_value_lollipop_plot.pdf", width = 8, height = 12)
ggplot(pvalue_filtered_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(xend = 0, yend = Description), size = 0.7) +  # Increase line size
  geom_point(aes(size = Count, color = pvalue)) +  # Increase dot size
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +  # Set minimum and maximum size for dots
  labs(title = "Enrichment Analysis Lollipop Plot",
       x = "Gene Count",
       y = "Enriched Pathway",
       color = "P-Value",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

# Interactive p-value lollipop plot
p <- ggplot(pvalue_filtered_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(xend = 0, yend = Description), size = 0.7) +  # Increase line size
  geom_point(aes(size = Count, color = pvalue)) +  # Increase dot size
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 7)) +  # Set minimum and maximum size for dots
  labs(title = "Enrichment Analysis Lollipop Plot",
       x = "Gene Count",
       y = "Enriched Pathway",
       color = "P-Value",
       size = "Gene Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = pvalue_text_size),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

# Convert the ggplot2 plot to an interactive plotly plot
interactive_plot <- ggplotly(p)

# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, "p_value_lollipop_plot.html")

# Network graph based on p-value
genes_pvalue <- unique(unlist(strsplit(pvalue_filtered_data$geneID, "/")))
pathways_pvalue <- unique(pvalue_filtered_data$Description)

# Create a data frame for edges (gene-pathway relationships)
edges_pvalue <- pvalue_filtered_data %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(from = geneID, to = Description)

# Create a data frame for nodes
nodes_pvalue <- data.frame(id = c(genes_pvalue, pathways_pvalue),
                           label = c(genes_pvalue, pathways_pvalue),
                           group = c(rep("gene", length(genes_pvalue)), rep("pathway", length(pathways_pvalue))))

# Create the visNetwork object for p-value filtered data
network_pvalue <- visNetwork(nodes_pvalue, edges_pvalue, width = "100%", height = "800px") %>%
  visNodes(shape = "dot", size = 10) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 123)

# Save the interactive plot as an HTML file
visSave(network_pvalue, file = "pvalue_cnetplot_interactive.html")

# Network graph based on adjusted p-value
genes_padjust <- unique(unlist(strsplit(padjust_filtered_data$geneID, "/")))
pathways_padjust <- unique(padjust_filtered_data$Description)

# Create a data frame for edges (gene-pathway relationships)
edges_padjust <- padjust_filtered_data %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(from = geneID, to = Description)

# Create a data frame for nodes
nodes_padjust <- data.frame(id = c(genes_padjust, pathways_padjust),
                            label = c(genes_padjust, pathways_padjust),
                            group = c(rep("gene", length(genes_padjust)), rep("pathway", length(pathways_padjust))))

# Create the visNetwork object for adjusted p-value filtered data
network_padjust <- visNetwork(nodes_padjust, edges_padjust, width = "100%", height = "800px") %>%
  visNodes(shape = "dot", size = 10) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 123)

# Save the interactive plot as an HTML file
visSave(network_padjust, file = "padjust_cnetplot_interactive.html")
