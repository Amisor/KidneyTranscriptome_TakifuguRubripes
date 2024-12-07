library(tidyverse)
library(topGO)

# ---------------------------Read files ---------------------------
path <- "/Users/ivanas.o/Desktop/MICB405_Docs/Final_Project/DESeq_topGO/"

# Upregulated and downregulated genes
up_gene <- read.csv(paste0(path,"results_filtered_up.csv"))
down_gene <- read.csv(paste0(path,"results_filtered_down.csv"))

# Load the mapping file 
geneID2GO <- readMappings(paste0(path, "takifugu_rubripes_GOIDs.tsv"))
UniGene <- names(geneID2GO)

# Extracting gene IDs 
upreg_genes <- as.character(up_gene$gene_id)
downreg_genes <- as.character(down_gene$gene_id)

# Factor names and set names 
up_gene_list <- factor(as.integer(UniGene %in% upreg_genes))
down_gene_list <- factor(as.integer(UniGene %in% downreg_genes))
names(up_gene_list) <- UniGene
names(down_gene_list) <- UniGene

# ---------------------------Go data ------------------------------------
# Build a GOdata object in topGO for upregulated
up_GO_data <- new("topGOdata", 
             description = "low", 
             ontology = "BP", 
             allGenes = up_gene_list,
             annot = annFUN.gene2GO,
             gene2GO = geneID2GO)

down_GO_data <- new("topGOdata",
               description = "low",
               ontology = "BP",
               allGenes = down_gene_list,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

# Run Fisher's test with the weight01 algorithm
up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")

down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

# Results summary for upregulated
up_GO <- GenTable(up_GO_data,
                     weight01 = up_result,
                     orderBy = "up_result",
                     ranksOf = "up_result",
                     topNodes = 50,
                     numChar = 1000)

# Results summary for down regulated
down_GO <- GenTable(down_GO_data,
                       weight01 = down_result,
                       orderBy = "down_result",
                       ranksOf = "down_result",
                       topNodes = 50,
                       numChar = 1000)

# -------------------- Visualizations Down -------------------
down_GO_filtered <- down_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

# First, let's arrange the data based on the enrichment ratio. 
down_GO_filtered_arranged <- down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Now let's extract the order of the term column
order_term <- down_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(title = "GO Enrichment Analysis for Downregulated Genes", x = "GO Term Description", y = "Enrichment Ratio", 
       color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1.10), breaks = seq(0, 1, 0.25), expand = c(0, 0.02)) # this changes the scale of the axes

#------------------- Visualizations Up ----------
up_GO_filtered <- up_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

# First, let's arrange the data based on the enrichment ratio. 
up_GO_filtered_arranged <- up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Now let's extract the order of the term column
order_term <- up_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(title = "GO Enrichment Analysis for Upregulated Genes", x = "GO Term Description", y = "Enrichment Ratio", 
       color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1.10), breaks = seq(0, 1, 0.25), expand = c(0, 0.02)) # this changes the scale of the axes


# -------------------Visualization Both -------------
# Add labels to upregulated and downregulated dataframes
up_GO <- up_GO %>%  mutate(up_down = "UP")

down_GO <- down_GO %>% mutate(up_down = "DOWN")

# Make a joined dataframe
joined_GO_filtered_arranged <- bind_rows(up_GO, down_GO) %>%
  filter(weight01 <= 0.05) %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term)) %>%
  head(n = 40)

# Extract the column order
order_term_joined <- joined_GO_filtered_arranged %>% 
  pull(Term)

# Plot joined dataframe
joined_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_point(aes(size= Significant)) +
  coord_flip() +
  scale_x_discrete(limits = order_term_joined) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(title = "GO Enrichment Analysis for DEGs", x = "GO Term Description", y = "Enrichment Ratio", 
       color = "P-value", size = "Number of Significant Genes") +
  theme(panel.border = element_rect(color = "black"), 
        panel.grid = element_line(colour = "grey96"), 
        strip.background = element_rect(colour = "black")) +
  scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.2), expand = c(0, 0.02)) +
  facet_grid(.~ up_down)
