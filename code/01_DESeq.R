library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# Function for reading tsv files
read_files <- function(path, file_tsv) {
  # path: directory where the file_tsv is
  # file_tsv: name of the tsv file - most match the name inside the directory
  columns_names <- c("gene_id", "total", "antisense", "sense")
  file <- read_tsv(paste0(path,file_tsv), col_names = columns_names, skip = 4)
  return(file)
}

path <- "/Users/ivanas.o/Desktop/MICB405_Docs/Final_Project/DESeq_topGO/"

# Low-temperature group (5°C)
Low1Reads <- read_files(path,"Low1Reads.tsv")
Low2Reads <- read_files(path, "Low2Reads.tsv")
Low3Reads <- read_files(path, "Low3Reads.tsv")

# Control group (23°C)
Control1Reads <- read_files(path,"Control1Reads.tsv")
Control2Reads <- read_files(path,"Control2Reads.tsv")
Control3Reads <- read_files(path,"Control3Reads.tsv")

# Create a dataframe with all the files
pufferfish_readcounts <- data.frame(row.names = Low1Reads$gene_id,
                                    Low1 = Low1Reads$sense,
                                    Low2 = Low2Reads$sense,
                                    Low3 = Low3Reads$sense,
                                    Control1 = Control1Reads$sense,
                                    Control2 = Control2Reads$sense,
                                    Control3 = Control3Reads$sense)

colnames(pufferfish_readcounts) # Order of columns is correct
class(pufferfish_readcounts) # It's a dataframe
dat_matrix<- as.matrix(pufferfish_readcounts) # Convert it to a matrix
head(dat_matrix) # Visualize the first counts
class(dat_matrix) # Class of the object

# Metadata - Pair each group
metadata <- data.frame(
  row.names = colnames(dat_matrix),
  condition = factor(c("LowTemp", "LowTemp", "LowTemp", "Control", "Control", "Control")),
  pair = factor(c("Pair1", "Pair2", "Pair3", "Pair1", "Pair2", "Pair3"))
)

head(metadata) # Visualize metadata

# Check if files are congruent
colnames(dat_matrix) == rownames(metadata)

# Check if converting values to integers will show any error
all_integers <- all(dat_matrix == as.integer(dat_matrix))

if (!all_integers) {
  print("There are non-integer values in the dat_matrix.")
} else {
  print("All values are integers.")
}

# We are ready to convert values to integers
dat_matrix <- as.matrix(dat_matrix)
mode(dat_matrix) <- "integer" 

# ----------------------------- DESeq2-------------------------------------
# Create DESeq2 object
dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix, #matrix 
                                     colData = metadata, #metadata file
                                     design = ~ pair + condition)

# Set control condition using the relevel function
dds_matrix$condition <- relevel(dds_matrix$condition, ref = "Control")

# Run DESeq2
dds <- DESeq(dds_matrix)

# Compare (contrast) LowTemp vs Control (obtain the fold-change)
results <- results(dds, contrast = c("condition", "LowTemp", "Control"))

# ----------------------------PCA --------------------------------------
# Log transform the count data and visualize it with PCA
rld <- rlog(dds)

# Extract PCA data
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)

# Add sample labels (use rownames from metadata or colnames of your dataset)
pca_data$Sample <- rownames(metadata)

# Plot PCA with non-overlapping labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = Sample)) +
  geom_point(size = 3) +  # Points for each sample
  geom_text_repel(size = 3, max.overlaps = 10) +  # Non-overlapping labels
  labs(
    title = "PCA of Gene Expression Data",
    x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
    y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")
  ) +
  theme_minimal()

#----------------------------Sample-to-sample distance ---------------------------- 
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, 
         main = "Sample-to-Sample Distances")

# ---------------------------Library size ------------------------------

# Total reads per sample
raw_reads <- colSums(counts(dds))
raw_reads

#Total reads across all samples
total_reads <- sum(raw_reads)
total_reads

# Group-specific reads
control_reads <- sum(raw_reads[grep("Control", names(raw_reads))])
low_temp_reads <- sum(raw_reads[grep("Low", names(raw_reads))])
control_reads
low_temp_reads

# ---------------------------Select results-----------------------------

# Names of the results that DESeq2 calculated 
# We don't need "pair_Pair2_vs_Pair1", "pair_Pair3_vs_Pair1"  
resultsNames(dds)

# Convert results to a dataframe to apply dplyr functions
# We only care about LowTemp_vs_Control
res <- results(dds, name = "condition_LowTemp_vs_Control") %>% as.data.frame()
res

# -------------------------- NA exploration --------------------------

# Number of genes with at least one NA value in any column
genes_with_na_any <- sum(apply(res, 1, function(x) any(is.na(x))))
print(paste("Number of genes with NA in any column:", genes_with_na_any))
nrow(res) # Number of genes
names(res) # Columns of DESeq results

# Number of genes with NA in the log2FoldChange column
genes_with_na_log2fc <- sum(is.na(res$log2FoldChange))
print(paste("Number of genes with NA in log2FoldChange:", genes_with_na_log2fc))

# Number of genes with NA in the pvalue column
genes_with_na_pvalue <- sum(is.na(res$"pvalue"))
print(paste("Number of genes with NA in pvalue:", genes_with_na_pvalue))

# Number of genes with NA in the padj column
genes_with_na_padj <- sum(is.na(res$"padj"))
print(paste("Number of genes with NA in padj:", genes_with_na_padj))

# Number of genes with NA in both log2FoldChange and pvalue
genes_with_na_log2fc_and_pvalue <- sum(is.na(res$log2FoldChange) & is.na(res$pvalue))
print(paste("Genes with NA in log2FoldChange and pvalue:", genes_with_na_log2fc_and_pvalue))

# Number of genes with NA in both log2FoldChange and padj
genes_with_na_log2fc_and_padj <- sum(is.na(res$log2FoldChange) & is.na(res$padj))
print(paste("Genes with NA in log2FoldChange and padj:", genes_with_na_log2fc_and_padj))

# Number of genes with NA in both pvalue and padj
genes_with_na_pvalue_and_padj <- sum(is.na(res$pvalue) & is.na(res$padj))
print(paste("Genes with NA in pvalue and padj:", genes_with_na_pvalue_and_padj))

# Number of genes with NA in log2FoldChange, pvalue, and padj
genes_with_na_all <- sum(is.na(res$log2FoldChange) & is.na(res$pvalue) & is.na(res$padj))
print(paste("Genes with NA in log2FoldChange, pvalue, and padj:", genes_with_na_all))

# Check if genes with zero in baseMean are the same as those with NA in pvalue
baseMean_zero_equals_pvalue_na <- all((res$baseMean == 0) == is.na(res$pvalue), na.rm = TRUE)
print(paste("Are genes with zero in baseMean the same as those with NA in pvalue?", baseMean_zero_equals_pvalue_na))

# Check if genes with zero in baseMean are the same as those with NA in log2FoldChange
baseMean_zero_equals_log2fc_na <- all((res$baseMean == 0) == is.na(res$log2FoldChange), na.rm = TRUE)
print(paste("Are genes with zero in baseMean the same as those with NA in log2FoldChange?", baseMean_zero_equals_log2fc_na))

# -----------------------------Remove NAN -------------------------------

# Getting rid of NA values
res_no_NA <- res %>% drop_na()
head(res_no_NA)
dim(res_no_NA) #look at how many rows you filtered out!
dim(res)

dim(res_no_NA)/dim(res)*100 # percentage of genes with no NA values = 63.7 %

# --------------------------------MA Plot----------------------------
green_color <- "#02c7b2" # Extracted green
red_color <- "#ed5817"   # Extracted red

ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  scale_x_log10() +
  labs(x = "log Mean of Normalized Counts", y = "Log2 Fold Change",
       title = "MA Plot of Differential Gene Expression") +
  theme_minimal() +
  # Separate horizontal lines with different colors
  geom_hline(yintercept = 1, color = green_color, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = -1, color = red_color, linetype = "dashed", size = 0.8) +
  scale_color_manual(values = c("#45473d", "steelblue"), 
                     labels = c("Not Significant", "Significant")) +
  theme(legend.position = "top")

#-------------------------- Filter: padj and log2FoldChange -----------

# Filtering for adjusted p-values <=0.05
res_filtered <- res_no_NA %>% filter(padj <= 0.05)

head(res_filtered)
dim(res_filtered) # look at how many rows you filtered out!

# Pull genes with more than 2x higher/lower expression
res_filtered_final <- res_filtered %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1)

# Function to write csv files
file_csv <- function(d, path, name){
  # Add the gene_id as a column to the dataframe
  new_d <-  cbind(gene_id = rownames(d), d)
  # No need to store rownames hthat contain the gene ids because are now stored as a column
  write.csv(new_d, paste0(path,name), row.names = FALSE)
  }
  
# Save DEGs 
file_csv(res_filtered_final, path, "results_filtered_final.csv")

# ----------------------- Upregulated and Downregulated genes ---------------

# Upregulated genes from the most upregulated to the least upregulated
upregulated_genes <- res_filtered_final %>% 
  arrange(desc(log2FoldChange)) %>% # desc() function to organize the column in descending order
  filter(log2FoldChange >= 1)

# Save upregulated genes
file_csv(upregulated_genes, path, "results_filtered_up.csv")

# Downregulated genes from the most downregulated to the least downregulated
downregulated_genes <- res_filtered_final %>% 
  arrange(log2FoldChange)  %>%
  filter(log2FoldChange <= -1)

# Save downregulated genes
file_csv(downregulated_genes, path, "results_filtered_down.csv")

# -------------------- Top 10 upregulated and downregulated genes -----------
# Top 10 upregulated genes
top10_genes <- upregulated_genes %>% head(n = 10)
file_csv(top10_genes, path,"top_10_DEGs.csv")
top10_genes

# Top 10 downregulated genes
bot10_genes <- downregulated_genes %>% head(n = 10)
file_csv(bot10_genes, path,"bot_10_DEGs.csv")

# --------------------------------- Heatmap z-scores -------------------
# Top 10 upregulated genes
top10_genes <- read.csv(paste0(path,"top_10_DEGs.csv"))

# Top 10 downregulated genes
bot10_genes <- read.csv(paste0(path,"bot_10_DEGs.csv"))

# Variance stabilizing transformation
vsd <- assay(vst(dds))
# Standarize the rows (genes) by subtracting the mean and dividing by the standard deviation
Z <- t(scale(t(vsd)))

# Top 10 upregulated genes z-scores
pheatmap(Z[top10_genes$gene_id,],
             main = "Heatmap of z-scores: Top 10 Upregulated Genes")

# Top 10 downregulated genes z-scores
pheatmap(Z[bot10_genes$gene_id,],
         main = "Heatmap of z-scores: Top 10 Downregulated Genes")
