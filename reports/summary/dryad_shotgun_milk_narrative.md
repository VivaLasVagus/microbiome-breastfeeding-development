# Dryad Milk Microbiome — Narrative Summary

# Overview

This narrative summarizes the ecological and diversity patterns observed in the Dryad shotgun-derived 
human milk microbiome dataset. The analysis focuses on community structure (ecotypes), richness and 
evenness, and longitudinal dynamics across lactation months. These results establish a clear ecological 
framework for understanding how the milk microbiome evolves over time and how individual-specific 
factors shape microbial balance. This work forms the foundation for downstream integration with infant 
gut 16S data and feeding-type metadata.

---

# Key Findings

- **Three distinct ecological types** were identified using Bray–Curtis distances, PCoA, hierarchical clustering, and k-means silhouette scoring.
- **Ecotype 1** represents a balanced, moderate-richness community and is the most stable across lactation.
- **Ecotype 2** is dominated by one or two genera (often Streptococcus or Corynebacterium), resulting in lower evenness.
- **Ecotype 3** is an individualized outlier profile driven by a single participant.
- **Richness increases with lactation month**, suggesting gradual diversification of the milk microbiome over time.
- **Evenness remains stable across lactation**, indicating that genus balance is shaped more by individual-specific factors than by time.
- **Participant identity explains variation in evenness**, but **not richness**, based on mixed-effects modeling.
- **Low biomass characteristics** of milk (low richness, high variability, sensitivity to contamination) are reflected in the dataset’s patterns.
- These ecological and diversity insights provide a strong foundation for linking milk microbiome structure to infant gut development in future analyses.

---

# Methods Summary

### **Data Source**
Shotgun metagenomic genus-level abundance data from the Dryad human milk microbiome dataset.

### **Preprocessing**
- Imported into R and converted into a phyloseq object.
- Collapsed to genus level
- Metadata cleaned and standardized

### **Ecological Typing**
- Bray–Curtis distance matrix computed from genus-level relative abundances.
- Principal Coordinates Analysis (PCoA) performed to visualize community structure.
- Hierarchical clustering and k-means clustering applied.
- Optimal number of clusters determined using silhouette scores.
- Three ecological types identified and assigned to samples.

### **Alpha Diversity**
- Richness calculated as observed genera (presence/absence).
- Shannon diversity computed using `vegan::diversity`.
- Pielou’s evenness derived from Shannon / log(richness).
- Diversity metrics merged with metadata for visualization and modeling.

### **Longitudinal Analysis**
- Richness and evenness plotted across lactation months.
- Participant trajectories visualized to assess stability and individual variation.

### **Mixed-Effects Modeling**
- Richness model: `Observed ~ Age (months) + (1 | Participant)`
- Evenness model: `Pielou ~ Age (months) + (1 | Participant)`
- Richness increased significantly with lactation month.
- Evenness showed no significant change over time.
- Participant random intercept collapsed for richness (singular fit), indicating no baseline differences.
- Participant identity explained variation in evenness.

---

# Narrative

The human milk microbiome in this dataset organizes into three distinct ecological types, 
each reflecting a different pattern of community structure. Ecotype 1 represents a balanced community 
with moderate richness and relatively high evenness, and it is the most common and stable pattern across 
participants. Ecotype 2 is characterized by dominance from one or two genera—often Streptococcus or 
Corynebacterium—resulting in lower evenness and greater variability across lactation months. 
Ecotype 3 is an outlier profile driven by a single participant, highlighting the individualized 
nature of low-biomass milk samples.

Across all ecological types, richness increases steadily with lactation month, 
suggesting that the milk microbiome becomes slightly more diverse over time, potentially due to infant 
oral seeding or stabilization of maternal skin-associated microbes. In contrast, 
evenness remains stable across lactation, indicating that the balance of genera within each sample 
is more strongly shaped by individual-specific factors than by time. 
Mixed-effects modeling supports this interpretation: participant identity explains variation in 
evenness but not richness, and ecological types capture meaningful differences in community structure 
that persist across the lactation period.

Together, these findings provide a clear ecological framework for understanding the milk microbiome 
and set the stage for downstream integration with infant gut 16S data. By characterizing both 
structural patterns (ecotypes) and diversity dynamics (richness and evenness), this analysis 
establishes a foundation for exploring how early-life microbial exposures may influence infant gut 
development.
