# 01_load_dryad.R
# Purpose: Load and inspect Dryad shotgun metagenome dataset

library(tidyverse)
BiocManager::install("phyloseq")

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

# QC
sample_names(ps_dryad)[1:10]   # peek at sample names
taxa_names(ps_dryad)[1:10]     # peek at taxa names

# Total abundance per sample
sample_sums <- data.frame(
  Sample = sample_names(ps_dryad),
  TotalAbundance = sample_sums(ps_dryad)
)

sample_sums %>% arrange(TotalAbundance)

# Sparsity
abund <- otu_table(ps_dryad)
sparsity <- sum(abund == 0) / length(abund)
sparsity

# Prevalence of each taxon
prev <- apply(abund, 1, function(x) sum(x > 0))
prev_df <- data.frame(
  Taxon = taxa_names(ps_dryad),
  Prevalence = prev
)

prev_df %>% arrange(Prevalence) %>% head(10)
prev_df %>% arrange(desc(Prevalence)) %>% head(10)

# Top taxa by mean abundance
mean_abund <- apply(abund, 1, mean)
top_taxa <- data.frame(
  Taxon = taxa_names(ps_dryad),
  MeanAbundance = mean_abund
) %>% arrange(desc(MeanAbundance))

top_taxa %>% head(10)

# Presence/absence matrix
pa_mat <- abund_mat > 0

# Observed richness per sample
richness <- colSums(pa_mat)

richness_df <- data.frame(
  Sample = names(richness),
  Observed = richness
)

richness_df %>% arrange(Observed)

# Removing zero-richness sample
ps_dryad <- prune_samples(richness > 0, ps_dryad)


library(ggplot2)

# Bray–Curtis distance matrix
bc_dist <- distance(ps_dryad, method = "bray")

pcoa_bc <- ordinate(ps_dryad, method = "PCoA", distance = bc_dist)

pcoa_df <- as.data.frame(pcoa_bc$vectors)
pcoa_df$Sample <- rownames(pcoa_df)

pcoa_df <- left_join(pcoa_df, meta, by = "Sample")

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = factor(Participant))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray–Curtis) — Colored by Participant",
       color = "Participant")

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = factor(`Age (months)`))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray–Curtis) — Colored by Lactation Month",
       color = "Age (months)")

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2)) +
  geom_point(color = "grey60", size = 2) +
  geom_point(data = pcoa_df %>% filter(Sample == "MN2892_VT1303"),
             color = "red", size = 4) +
  theme_minimal() +
  labs(title = "PCoA — Highlighting High-Richness Outlier")

# Removing outlier
ps_dryad <- prune_samples(sample_names(ps_dryad) != "MN2892_VT1303", ps_dryad)

bc_dist <- distance(ps_dryad, method = "bray")

meta_ord <- meta[match(rownames(as.matrix(bc_dist)), meta$Sample), ]

library(vegan)

adonis2(
  bc_dist ~ Participant + `Age (months)`,
  data = meta_ord,
  permutations = 999
)

adonis2(
  bc_dist ~ Participant * `Age (months)`,
  data = meta_ord,
  permutations = 999
)

# Prevalence threshold: taxa present in at least 3 samples
prev <- apply(abund_mat, 1, function(x) sum(x > 0))
keep_taxa <- names(prev[prev >= 3])

ps_filt <- prune_taxa(keep_taxa, ps_dryad)
ps_filt

mean_abund <- apply(otu_table(ps_filt), 1, mean)
keep_taxa2 <- names(mean_abund[mean_abund >= 0.001])

ps_filt <- prune_taxa(keep_taxa2, ps_filt)
ps_filt

# Recheck sparsity
abund_filt <- otu_table(ps_filt)
sparsity_filt <- sum(abund_filt == 0) / length(abund_filt)
sparsity_filt

# Rerun ordination
bc_dist_filt <- distance(ps_filt, method = "bray")
pcoa_filt <- ordinate(ps_filt, method = "PCoA", distance = bc_dist_filt)

pcoa_df_filt <- as.data.frame(pcoa_filt$vectors)
pcoa_df_filt$Sample <- rownames(pcoa_df_filt)
pcoa_df_filt <- left_join(pcoa_df_filt, meta, by = "Sample")

ggplot(pcoa_df_filt, aes(x = Axis.1, y = Axis.2, color = factor(Participant))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Filtered PCoA (Bray–Curtis) — Colored by Participant")

# Rerun PERMANOVA
meta_ord_filt <- meta[match(rownames(as.matrix(bc_dist_filt)), meta$Sample), ]

adonis2(
  bc_dist_filt ~ Participant + `Age (months)`,
  data = meta_ord_filt,
  permutations = 999
)

# Collapse to genus level
ps_genus <- tax_glom(ps_filt, taxrank = "Genus")
ps_genus

# Top 10 genera
genus_abund <- apply(otu_table(ps_genus), 1, mean)
sort(genus_abund, decreasing = TRUE)[1:10]

# Melt genus table
genus_long <- otu_table(ps_genus) %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance")

# Add metadata
genus_long <- left_join(genus_long, meta, by = "Sample")

# Compute mean abundance per participant × genus
genus_participant <- genus_long %>%
  group_by(Participant, Genus) %>%
  summarize(MeanAbundance = mean(Abundance), .groups = "drop")

# Apply abundance threshold (e.g., ≥ 0.01)
genus_participant_thresh <- genus_participant %>%
  filter(MeanAbundance >= 0.01)

# Select top 5 genera per participant
top5_thresh <- genus_participant_thresh %>%
  group_by(Participant) %>%
  slice_max(order_by = MeanAbundance, n = 5)

# Plot
ggplot(top5_thresh,
       aes(x = reorder(Genus, MeanAbundance), y = MeanAbundance, fill = Genus)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Participant, scales = "free_y") +
  theme_minimal() +
  labs(title = "Top Genera per Participant (Threshold Applied)",
       x = "Genus",
       y = "Mean Relative Abundance")

genus_long <- otu_table(ps_genus) %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance")

genus_long <- left_join(genus_long, meta, by = "Sample")

# Identify top genera overall
top_genera <- genus_long %>%
  group_by(Genus) %>%
  summarize(MeanAbundance = mean(Abundance)) %>%
  arrange(desc(MeanAbundance)) %>%
  slice(1:5) %>%
  pull(Genus)

# Filter to top genera
genus_long_top <- genus_long %>%
  filter(Genus %in% top_genera)

# Plot
ggplot(genus_long_top,
       aes(x = `Age (months)`, y = Abundance, color = Genus)) +
  geom_line(aes(group = Sample), alpha = 0.4) +
  geom_smooth(se = FALSE, size = 1.2) +
  theme_minimal() +
  labs(title = "Longitudinal Trajectories of Top Genera Across Lactation",
       x = "Lactation Month",
       y = "Relative Abundance")

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
       aes(x = `Age (months)`,
           y = Abundance,
           color = Genus,
           group = interaction(Participant, Genus))) +
  geom_line(alpha = 0.7, size = 1.1) +
  geom_point(size = 2.5, alpha = 0.9) +
  facet_wrap(~ Participant, scales = "free_y") +
  theme_minimal() +
  labs(title = "Longitudinal Trajectories of Top 3 Genera per Participant",
       x = "Lactation Month",
       y = "Relative Abundance")


# Checking data for participants
table(meta$Participant)


# Comparing phyloseq and metadata
sample_names(ps_genus)
meta$Sample

setdiff(meta$Sample, sample_names(ps_genus))
setdiff(sample_names(ps_genus), meta$Sample)

meta %>% filter(Sample %in% c("MN2892_VT1303", "MN2336_VT1303"))

meta <- meta %>% 
  filter(Sample %in% sample_names(ps_genus))

# Build genus x sample matrix
genus_mat <- otu_table(ps_genus) %>%
  as.data.frame()

# Make sure rows = samples, columns = genera
genus_mat <- t(genus_mat)

# Adding participant and month metadata
meta_heat <- meta %>%
  select(Sample, Participant, `Age (months)`)

# Reorder metadata to match matrix row order
meta_heat <- meta_heat[match(rownames(genus_mat), meta_heat$Sample), ]

# Selecting top genera
top_genera <- colMeans(genus_mat) %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%        # adjust to 10, 15, 20 depending on clarity
  names()

genus_mat_top <- genus_mat[, top_genera]

# Remove genera with all zeros or all NA
genus_mat_top <- genus_mat_top[, colSums(genus_mat_top, na.rm = TRUE) > 0]

# Remove genera with zero variance across samples
genus_mat_top <- genus_mat_top[, apply(genus_mat_top, 2, var, na.rm = TRUE) > 0]

# Building heatmap
library(pheatmap)

pheatmap(genus_mat_top,
         scale = "row",                     # normalize per genus
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_row = meta_heat[, c("Participant", "Age (months)")],
         fontsize_row = 6,
         fontsize_col = 8,
         main = "Heatmap of Top Genera Across Samples")
# Build genus x sample matrix
genus_mat <- otu_table(ps_genus) %>%
  as.data.frame()

# Make sure rows = samples, columns = genera
genus_mat <- t(genus_mat)

# Adding participant and month metadata
meta_heat <- meta %>%
  select(Sample, Participant, `Age (months)`)

# Reorder metadata to match matrix row order
meta_heat <- meta_heat[match(rownames(genus_mat), meta_heat$Sample), ]

# Selecting top genera
top_genera <- colMeans(genus_mat) %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%        # adjust to 10, 15, 20 depending on clarity
  names()

genus_mat_top <- genus_mat[, top_genera]

# Remove genera with all zeros or all NA
genus_mat_top <- genus_mat_top[, colSums(genus_mat_top, na.rm = TRUE) > 0]

# Remove genera with zero variance across samples
genus_mat_top <- genus_mat_top[, apply(genus_mat_top, 2, var, na.rm = TRUE) > 0]

library(tidyverse)

heat_long <- genus_mat_top %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Genus", values_to = "Abundance")

ggplot(heat_long, aes(x = Genus, y = Sample, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Heatmap of Top Genera Across Samples",
       x = "Genus",
       y = "Sample")

