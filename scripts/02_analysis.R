###############################################
# 02_analysis.R
# Purpose: Ecological clustering, ordination overlays,
#          differential abundance, participant microbiome types
###############################################

###############################################
# Reproducibility
###############################################
set.seed(2026)
sessionInfo()

###############################################
# Load libraries
###############################################
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(cluster)
library(factoextra)

###############################################
# Load processed phyloseq object
###############################################
ps_genus <- readRDS("data/processed/ps_genus.rds")

meta <- data.frame(sample_data(ps_genus))
meta$Participant <- factor(meta$Participant)

meta <- data.frame(sample_data(ps_genus)) %>%
  rename(`Age (months)` = Age..months.)


###############################################
# 1. Bray–Curtis distance + PCoA overlays
###############################################

bc_dist <- distance(ps_genus, method = "bray")
pcoa <- ordinate(ps_genus, method = "PCoA", distance = bc_dist)

pcoa_df <- as.data.frame(pcoa$vectors) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by = "Sample")

###############################################
# PCoA colored by Participant
###############################################

p_pcoa_participant <- ggplot(pcoa_df,
                             aes(Axis.1, Axis.2, color = Participant)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_discrete() +
  labs(title = "PCoA (Bray–Curtis) — Participant",
       x = "Axis 1", y = "Axis 2")

ggsave("results/figures/pcoa_participant_clean.png",
       p_pcoa_participant, width = 8, height = 6, dpi = 300)

###############################################
# PCoA colored by Lactation Month
###############################################

p_pcoa_month <- ggplot(pcoa_df,
                       aes(Axis.1, Axis.2, color = `Age (months)`)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_viridis_c() +
  labs(title = "PCoA (Bray–Curtis) — Lactation Month",
       x = "Axis 1", y = "Axis 2")

ggsave("results/figures/pcoa_month_clean.png",
       p_pcoa_month, width = 8, height = 6, dpi = 300)

###############################################
# 2. Ecological clustering (participant-level)
###############################################

# Build participant × genus matrix
genus_mat <- as.matrix(t(otu_table(ps_genus)))

# Compute participant-level mean abundance
genus_participant <- genus_mat %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by = "Sample") %>%
  group_by(Participant) %>%
  summarize(across(where(is.numeric), mean), .groups = "drop")

# Remove Participant column for clustering
clust_mat <- genus_participant %>%
  select(-Participant, -Age (months)) %>%   # remove metadata
  as.matrix()

rownames(clust_mat) <- genus_participant$Participant


###############################################
# Hierarchical clustering
###############################################

hc <- hclust(dist(clust_mat), method = "ward.D2")

p_hclust <- fviz_dend(hc,
                      k = 4,                      # number of ecological types
                      cex = 1.2,
                      color_labels_by_k = TRUE,
                      rect = TRUE,
                      main = "Participant Ecological Clusters")

ggsave("results/figures/ecological_clusters.png",
       p_hclust, width = 8, height = 6, dpi = 300)

###############################################
# 3. K-means clustering + silhouette analysis
###############################################

sil_scores <- data.frame()

for (k in 2:8) {
  km <- kmeans(clust_mat, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(clust_mat))
  sil_scores <- rbind(sil_scores,
                      data.frame(k = k,
                                 mean_sil = mean(ss[, 3])))
}

write_csv(sil_scores, "results/tables/silhouette_scores.csv")

p_sil <- ggplot(sil_scores, aes(k, mean_sil)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Silhouette Scores for K-means Clustering",
       x = "Number of clusters (k)",
       y = "Mean silhouette width")

ggsave("results/figures/silhouette_scores.png",
       p_sil, width = 8, height = 6, dpi = 300)

###############################################
# 4. Differential abundance across lactation
###############################################

genus_long <- genus_mat %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Genus", values_to = "Abundance") %>%
  left_join(meta, by = "Sample")

# Simple linear model per genus
da_results <- genus_long %>%
  group_by(Genus) %>%
  do({
    fit <- lm(Abundance ~ `Age (months)`, data = .)
    data.frame(
      Genus = unique(.$Genus),
      p_value = summary(fit)$coefficients[2, 4],
      slope = summary(fit)$coefficients[2, 1]
    )
  })

da_results <- da_results %>%
  arrange(p_value)

write_csv(da_results, "results/tables/differential_abundance_lactation.csv")

###############################################
# 5. Participant-level heatmap (ecological types)
###############################################

library(pheatmap)

pheatmap(clust_mat,
         scale = "row",
         clustering_method = "ward.D2",
         main = "Participant × Genus Heatmap",
        filename = "results/figures/participant_genus_heatmap.png",
        width = 10,
        height = 12)

###############################################
# End of 02_analysis.R
###############################################
