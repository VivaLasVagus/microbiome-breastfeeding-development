###############################################
# Script 01 — MIMBES Genus-Level Import
###############################################

library(tidyverse)
library(phyloseq)
library(janitor)

# Load genus abundance data
load("data/raw/mimbe_milk_16s/genusdata.Rdata")  
# This usually loads an object named something like 'genusdata'

# Load metadata
meta <- read.csv("data/raw/mimbe_milk_16s/metadata.csv") %>% clean_names()

# Inspect objects
str(genusdata)
head(meta)

