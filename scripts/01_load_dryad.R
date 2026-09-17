# 01_load_dryad.R
# Purpose: Load and inspect Dryad shotgun metagenome dataset

library(tidyverse)


# Paths
taxa_path <- "data/raw/dryad_milk_shotgun/Milk_Metagenome_Taxa_Table.txt"
meta_path <- "data/raw/dryad_milk_shotgun/Sample_Age_Participant.txt"

# Load data
taxa <- read.delim(taxa_path, check.names = FALSE)
meta <- read.delim(meta_path, check.names = FALSE)

# Peek at structure
glimpse(taxa)
glimpse(meta)

# Separate taxonomy and abundance
tax_table_df <- taxa[, 1:8]
abund_df <- taxa[, 9:(ncol(taxa)-2)]

# Convert abundance to matrix
abund_mat <- as.matrix(abund_df)
rownames(abund_mat) <- taxa$Name

# Convert taxonomy to matrix
tax_mat <- as.matrix(tax_table_df)
rownames(tax_mat) <- taxa$Name

# Align metadata rownames
meta2 <- meta
rownames(meta2) <- meta2$Sample

# Build phyloseq object
library(phyloseq)

ps_dryad <- phyloseq(
  otu_table(abund_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta2)
)

ps_dryad
