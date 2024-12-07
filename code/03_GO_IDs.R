library(tidyverse)

path <- "/Users/ivanas.o/Desktop/MICB405_Docs/Final_Project/DESeq_topGO/"

# Download data from NCBI
gaf_data <- read.table(paste0(path,"GCF_901000725.2_fTakRub1.2_gene_ontology.gaf"), 
                       sep = "\t", header = FALSE, 
                       comment.char = "!", quote = "", fill = TRUE)

# Set columns names manually
colnames(gaf_data) <- c("DB", "GeneID", "Symbol", "Qualifier", "GO_ID", "Reference", 
                        "Evidence_Code", "With_From", "Aspect", "Gene_Name", "Gene_Synonym", 
                        "Type", "Taxon", "Date", "Assigned_By", "Annot_Ext", 
                        "Gene_Product_Form_ID")

# Focus on the Symbol - corresponding to gene_id of the ref genome and the GO_ID
gaf_data_f <- gaf_data[,c("Symbol", "GO_ID")]

# Create a collapse dataframe. Have one gene_id and all its GO_IDs
go_terms <- gaf_data_f %>% 
  group_by(Symbol) %>% 
  summarize(GO_ID = paste(GO_ID, collapse = ", "))

# Rename columns
names(go_terms) <- c("GeneID", "GOTerms")

# Remove undesired characters from the dataframe columns, like leading/trailing quotes
go_terms$GeneID <- gsub('(^"|"$)', '', go_terms$GeneID) 
go_terms$GOTerms <- gsub('(^"|"$)', '', go_terms$GOTerms) 

# Write the gene ids and their associated go terms
write.table(go_terms, file = paste0(path, "/takifugu_rubripes_GOIDs.tsv"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

