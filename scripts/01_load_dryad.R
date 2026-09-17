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

