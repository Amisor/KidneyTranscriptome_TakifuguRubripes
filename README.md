# KidneyTranscriptome_TakifuguRubripes

This repository contains the code and analyses for the project **Low-temperature Stress Alters Gene Expression of the Kidney Transcriptome in Takifugu rubripes**.

## Contents

### Analyses
- FastQC reports: Quality control reports for the raw sequencing data.
- DESeq2 analysis: Identification of differentially expressed genes (DEGs).
- TopGO enrichment analysis: Identification of key biological processes associated with the DEGs.

### Directories

#### `data`
Contains:
- Read counts generated from STAR:  
  - `Low1Reads.tsv`, `Low2Reads.tsv`, `Low3Reads.tsv`  
  - `Control1Reads.tsv`, `Control2Reads.tsv`, `Control3Reads.tsv`
- Gene ontology annotation files:  
  - `.gaf` file: `GCF_901000725.2_fTakRub1.2_gene_ontology.gaf`  
  - `.tsv` file linking gene IDs to GO terms: `takifugu_rubripes_GOIDs.tsv`
- Differentially expressed genes (DEGs):  
  - All DEGs: `results_filtered_final.csv`  
  - Upregulated genes: `results_filtered_up.csv`  
  - Downregulated genes: `results_filtered_down.csv`
- Top 10 DEGs:  
  - Top 10 upregulated: `top_10_DEGs.csv`  
  - Top 10 downregulated: `bot_10_DEGs.csv`

#### `code`
Includes:
- `01_DESeq.R`: DESeq2 analysis for identifying DEGs, PCA, sample-to-sample distances, MA plot, and heatmap of z-scores for the top 10 DEGs.
- `02_DEGs.R`: Visualization of overall DEGs divided into control and low-temperature groups, upregulated and downregulated genes, and plots for the top 10 upregulated and downregulated genes.
- `03_GO_IDs.R`: Script for generating the GO terms mapping file.
- `04_topGO.R`: GO enrichment analysis using TopGO.

#### `FastQC`
Contains:
- PDF and HTML versions of FastQC reports for the different samples.
