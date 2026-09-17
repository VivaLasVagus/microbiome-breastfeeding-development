###############################################
# 01_load_dryad_clean.R
# Purpose: Load, QC, filter, analyze Dryad shotgun metagenome dataset
###############################################

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

###############################################
# 1. Load data
###############################################

taxa_path <- "data/raw/dryad_milk_shotgun/Milk_Metagenome_Taxa_Table.txt"
meta_path <- "data/raw/dryad_milk_shotgun/Sample_Age_Participant.txt"

taxa <- read.delim(taxa_path, check.names = FALSE)
meta <- read.delim(meta_path, check.names = FALSE)

###############################################
# 2. Prepare taxonomy + abundance matrices
###############################################

tax_table_df <- taxa[, 1:8]
abund_df <- taxa[, 9:(ncol(taxa)-2)]

abund_mat <- as.matrix(abund_df)
rownames(abund_mat) <- taxa$Name

tax_mat <- as.matrix(tax_table_df)
rownames(tax_mat) <- taxa$Name

###############################################
# 3. Sync metadata to abundance table
###############################################

meta <- meta %>% filter(Sample %in% colnames(abund_mat))
rownames(meta) <- meta$Sample

###############################################
# 4. Build phyloseq object
###############################################

ps_dryad <- phyloseq(
  otu_table(abund_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta)
)

###############################################
# 5. QC: richness, sparsity, prevalence
###############################################

pa_mat <- abund_mat > 0
richness <- colSums(pa_mat)

ps_dryad <- prune_samples(richness > 0, ps_dryad)

###############################################
# 6. Ordination + PERMANOVA
###############################################

bc_dist <- distance(ps_dryad, method = "bray")
pcoa_bc <- ordinate(ps_dryad, method = "PCoA", distance = bc_dist)

pcoa_df <- as.data.frame(pcoa_bc$vectors) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by = "Sample")

ggplot(pcoa_df, aes(Axis.1, Axis.2, color = Participant)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray–Curtis) — Participant")

###############################################
# 7. Prevalence + abundance filtering
###############################################

abund <- otu_table(ps_dryad)
prev <- apply(abund, 1, function(x) sum(x > 0))
keep_taxa <- names(prev[prev >= 3])

ps_filt <- prune_taxa(keep_taxa, ps_dryad)

mean_abund <- apply(otu_table(ps_filt), 1, mean)
keep_taxa2 <- names(mean_abund[mean_abund >= 0.001])

ps_filt <- prune_taxa(keep_taxa2, ps_filt)

###############################################
# 8. Collapse to genus level
###############################################

ps_genus <- tax_glom(ps_filt, taxrank = "Genus")

###############################################
# 9. Longitudinal top 3 genera per participant
###############################################

meta$Participant <- as.factor(meta$Participant)

genus_long <- otu_table(ps_genus) %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance") %>%
  left_join(meta, by = "Sample")

top3 <- genus_long %>%
  group_by(Participant, Genus) %>%
  summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  group_by(Participant) %>%
  slice_max(order_by = MeanAbundance, n = 3)

genus_long_top3 <- genus_long %>%
  semi_join(top3, by = c("Participant", "Genus"))

ggplot(genus_long_top3,
       aes(`Age (months)`, Abundance, color = Genus,
           group = interaction(Participant, Genus))) +
  geom_line(alpha = 0.7, size = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(~ Participant, scales = "free_y") +
  theme_minimal() +
  labs(title = "Longitudinal Trajectories of Top 3 Genera per Participant")

###############################################
# 10. Heatmap
###############################################

genus_mat <- as.matrix(t(otu_table(ps_genus)))

meta_heat <- meta %>%
  select(Sample, Participant, `Age (months)`) %>%
  filter(Sample %in% rownames(genus_mat)) %>%
  slice(match(rownames(genus_mat), Sample))

top_genera <- colMeans(genus_mat) %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%
  names()

genus_mat_top <- genus_mat[, top_genera]

heat_long <- genus_mat_top %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Genus", values_to = "Abundance")

ggplot(heat_long, aes(Genus, Sample, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Heatmap of Top Genera Across Samples")
