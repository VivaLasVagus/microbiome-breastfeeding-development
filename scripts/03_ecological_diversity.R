###############################################
# 03_ecological_diversity.R
# Purpose: analyzing evenness and richness,
#          ecological-type comparison
###############################################

###############################################
# Reproducibility
###############################################
set.seed(2026)
sessionInfo()

###############################################
# Load libraries
###############################################
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)

####################################################
# Load processed phyloseq object & create data frame
####################################################
ps_genus <- readRDS("data/processed/ps_genus.rds")

meta <- data.frame(sample_data(ps_genus)) %>%
  rename(`Age (months)` = Age..months.)

meta$Participant <- factor(meta$Participant)

###############################################
# Compute richness & evenness
###############################################

abund_mat <- as.data.frame(t(otu_table(ps_genus)))

richness <- rowSums(abund_mat > 0)
shannon <- vegan::diversity(abund_mat, index = "shannon")
evenness <- shannon / log(richness)
evenness[richness == 0] <- NA

alpha_df <- data.frame(
  SampleID = rownames(abund_mat),
  Observed = richness,
  Shannon = shannon,
  Pielou = evenness
)

# Add metadata
alpha_df <- cbind(meta, alpha_df)

###############################################
# Visualize richness & evenness over time
###############################################

# Richness boxplot
p_rich_box <- ggplot(alpha_df, aes(x = `Age (months)`, y = Observed)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "Richness Across Lactation Months")

ggsave("results/figures/richness_across_months.png",
       p_rich_box, width = 8, height = 6, dpi = 300)

# Evenness boxplot
p_even_box <- ggplot(alpha_df, aes(x = `Age (months)`, y = Pielou)) +
  geom_boxplot(fill = "tan") +
  theme_minimal() +
  labs(title = "Evenness Across Lactation Months")

ggsave("results/figures/evenness_across_months.png",
       p_even_box, width = 8, height = 6, dpi = 300)

# Participant richness trajectories
p_rich_traj <- ggplot(alpha_df, aes(x = `Age (months)`, y = Observed,
                                    group = Participant, color = Participant)) +
  geom_line(alpha = 0.6) +
  geom_point() +
  theme_minimal() +
  labs(title = "Participant Richness Trajectories")

ggsave("results/figures/richness_trajectories.png",
       p_rich_traj, width = 8, height = 6, dpi = 300)


###############################################
# Mixed-effect statistical testing
###############################################

# Richness model
rich_mod <- lmer(Observed ~ `Age (months)` + (1 | Participant), data = alpha_df)
summary(rich_mod)

# Evenness model
even_mod <- lmer(Pielou ~ `Age (months)` + (1 | Participant), data = alpha_df)
summary(even_mod)

