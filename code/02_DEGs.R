library(tidyverse)
library(ggplot2)
library(patchwork)

# --------------------------------- Read files --------------------------------
path <- "/Users/ivanas.o/Desktop/MICB405_Docs/Final_Project/DESeq_topGO/"

# Function for reading tsv files
read_files <- function(path, file_csv) {
  # path: directory where the file_csv is
  # file_csv: name of the tsv file - most match the name inside the directory
  columns_names <- c("gene_id", "total", "antisense", "sense")
  file <- read.csv(paste0(path,file_csv))
  return(file)
}

# DEGs
results_filtered <- read_files(path,"results_filtered_final.csv")

# Upregulated genes
up <- read_files(path,"results_filtered_up.csv")

# Downregulated genes
down <- read_files(path,"results_filtered_down.csv")

# Top 10 upregulated genes
top10_genes <- read_files(path,"top_10_DEGs.csv")

# Top 10 downregulated genes
bot10_genes <- read_files(path,"bot_10_DEGs.csv")

# Normalized counts 
normalized_counts <- read_files(path,"normalized_counts.csv")
normalized_counts$condition <- factor(normalized_counts$condition, levels = c("LowTemp", "Control")) 

# --------------------------Plot general results -----------------
genes_tot <- length(results_filtered$gene_id) # Total DEGs
genes_up <- length(up$gene_id) # Number of upregulated genes
genes_down <- length(down$gene_id) # Number of upregulated genes

categories <- c("Up", "Down", "All")
values <- c(genes_up, genes_down, genes_tot)

# Create the barplot
barplot(values,
        names.arg = categories,
        col = c(green_color, red_color, "gray"), # Set colors for Up, Down, and All
        ylim = c(0, max(values) + 100), # Adding some space above the tallest bar
        ylab = "Number of DEGs",
        xlab = "Low Temperature (5°C) vs Control (23°C)",
        main = "Differentially Expressed Genes (DEGs)",
        font.main = 1,
        border = "black",
        las = 1, # Rotate axis labels for better visibility
        cex.names = 1.2) # Adjust label size for better visibility

# Add numbers above the bars
text(x = 1:length(values),
     y = values,
     labels = values,
     pos = 3, # Position text above the bar
     cex = 1) # Text size

# ------------------------ Top 10 upregulated and downregulated genes ---------------------

# Combine top and bottom genes into a single data frame
combined_genes <- rbind(top10_genes, bot10_genes)
combined_genes$type <- ifelse(combined_genes$log2FoldChange > 0, 
                              "Upregulated", "Downregulated")

ggplot(combined_genes, aes(x = reorder(gene_id, log2FoldChange), 
                           y = log2FoldChange, color = type, size = -log10(padj))) +
  geom_point() +
  coord_flip() +
  labs(x = "Gene ID", y = "Log2 Fold Change", 
       title = "Top 10 Upregulated and Downregulated Genes") +
  theme_minimal()

# ----------------------------- Function to plot nomalized counts -----------

norm_counts_plot <- function(normalized_counts, top_data){
  # Store individual plots in a list
  plot_list <- list()
  
  for (i in seq_along(top_data$gene_id)) {
    # Plot per sample
    gene_of_interest <- top_data$gene_id[i] # Indenitfy gene
    plot <- normalized_counts %>% # Generate the plot
      filter(gene_id == gene_of_interest) %>%
      ggplot(aes(x = Sample, y = count, color = Sample)) + 
      geom_jitter(width = 0.2, size = 3) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1), # rotate lables
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +  # Add extra space above and below
      ggtitle(gene_of_interest)
    
    # Add plot to the list
    plot_list[[i]] <- plot
  }
  return(plot_list)
}

# -----------------------------Plot top 10 upregulated genes normalized counts --------------
# Get a list of all the plots of the top 10 upregulated genes
up_plot_list <- norm_counts_plot(normalized_counts, top10_genes)

# Combine all plots into a single figure (2 columns, 5 rows) for upregulated genes
up_combined_plot <- wrap_plots(up_plot_list, ncol = 2, guides = "collect") +
  plot_layout(guides = "collect") +  # Ensure shared legend
  plot_annotation(title = "Top 10 Upregulated Genes: Normalized Counts per Sample") 

# Print plot for upregulated genes
print(up_combined_plot)

ggsave(
  filename = paste0(path, "up_plot.png"),
  plot = up_combined_plot,
  width = 10,    # Width of the image (in inches)
  height = 12,   # Height of the image (in inches)
  dpi = 300      # Resolution in dots per inch (300 for high quality)
)

# -----------------------------Plot top 10 downregulated genes normalized counts --------------
# Get a list of all the plots of the top 10 downregulated genes
down_plot_list <- norm_counts_plot(normalized_counts, bot10_genes)

# Combine all plots into a single figure (2 columns, 5 rows) for downregulated genes
down_combined_plot <- wrap_plots(down_plot_list, ncol = 2, guides = "collect") +
  plot_layout(guides = "collect") +  # Ensure shared legend
  plot_annotation(title = "Top 10 Downregulated Genes: Normalized Counts per Sample") 

# Print plot for upregulated genes
print(down_combined_plot)

ggsave(
  filename = paste0(path, "down_plot.png"),
  plot = down_combined_plot,
  width = 10,    # Width of the image (in inches)
  height = 12,   # Height of the image (in inches)
  dpi = 300      # Resolution in dots per inch (300 for high quality)
)
