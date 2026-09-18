# Ecological Types — Summary of Findings

This document summarizes the ecological types identified in the human milk microbiome dataset and key diversity results from Script 03.

---

## Overview

Using genus-level abundance data from the Dryad milk microbiome dataset, three ecological types were identified through Bray–Curtis PCoA, hierarchical clustering, and k-means silhouette scoring. These ecological types represent distinct community structures within human milk samples.

---

## Ecological Types

### **Ecotype 1 — Balanced Community**
- Moderate richness (median ~5 genera)
- High evenness (Pielou ~0.6–0.8)
- No single genus dominates
- Represents the majority of participants
- Likely reflects typical lactation-associated microbial exposure

### **Ecotype 2 — Streptococcus/Corynebacterium-Dominant**
- Lower evenness (skewed distribution)
- One or two genera dominate the community
- May reflect maternal skin/oral transfer patterns
- Shows more variability across lactation months

### **Ecotype 3 — Outlier Profile**
- Driven by a single participant
- Unique genus composition
- Important for understanding individual variation and potential contamination or biological uniqueness

---

## Diversity Results

### **Richness (Observed Genera)**
- Milk samples show low richness (2–12 genera)
- Median richness ~5 genera
- Consistent with low-biomass nature of human milk

### **Evenness (Pielou’s Index)**
- Evenness ranges from ~0.2 to 1.0
- Median ~0.65
- High evenness is typical in low-biomass samples where few genera share similar proportions

### **Longitudinal Patterns**
- Some participants show stable richness/evenness across lactation months
- Others show variability, likely due to sampling differences, maternal factors, or low-biomass sensitivity

---

## Ecotype Comparisons

Mixed‑effects models indicate:

Ecotypes differ in richness and evenness  
(This comes from your boxplots + clustering results, not the mixed models.)

Participant identity explains variation in evenness,
but not richness  
(Richness model had a singular fit → random intercept collapsed to zero.)

Richness increases with lactation month,
regardless of ecotype
(Slope = +0.56 genera per month, p = 0.002)

Evenness does not change with lactation month  
(Slope ~0, p = 0.565)

Ecotype 2 tends to have lower evenness  
(Matches your ecological‑type boxplots.)

Ecotype 1 appears most stable across time  
(Matches your longitudinal trajectories.)

---

## Next Steps

- Begin infant gut 16S MIMBES pipeline  
- Test correlations between milk ecotypes and infant gut profiles  
- Integrate feeding-type metadata  
- Build models linking milk microbiome → infant gut → developmental outcomes  

