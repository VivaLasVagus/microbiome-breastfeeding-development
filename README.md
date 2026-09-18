# microbiome-breastfeeding-development

A reproducible analysis exploring how breastfeeding shapes early-life microbiome ecology — beginning with the human milk microbiome — and how these microbial patterns may relate to infant gut development. This project is part of my learning and portfolio-building journey in microbiome science, maternal/child health, and early-life development.

---

## Project Goals

- Characterize the human milk microbiome using shotgun-derived genus-level data  
- Identify ecological types (community structure patterns) in milk samples  
- Examine richness, evenness, and longitudinal stability across lactation months  
- Build reproducible R workflows for microbiome analysis (phyloseq, vegan, tidyverse)  
- Prepare for downstream integration with infant gut 16S data (MIMBES pipeline)  
- Develop statistical models linking microbial features to feeding type and developmental outcomes  
- Practice reproducible workflows using GitHub, R, Quarto, and tidy data principles  

---

## Repository Structure

data/
raw/            # Unmodified source data
processed/      # Cleaned and transformed phyloseq objects

scripts/
01_qc_preprocessing.R        # Data ingestion, filtering, genus-level collapse
02_ecological_typing.R       # Bray–Curtis PCoA, clustering, silhouette scores
03_ecological_diversity.R    # Richness, evenness, stability, ecotype comparisons

Future: 
04_16S_MIMBES_pipeline.R
05_visualizations.R
06_stats_models.R

notebooks/
exploratory_analysis.qmd     # Quarto notebooks for exploratory work
functional_pathways.qmd

reports/
microbiome_summary.md        # Narrative summaries
developmental_outcomes.md
summary/
ecological_types.md        # Ecotype definitions & results summary

dashboards/
powerbi/                     # Power BI files for interactive visualization

results/                  #Plots and .csv files from R


---

## Tools & Technologies

- **R** (tidyverse, phyloseq, vegan, microbiome, ggplot2)  
- **Quarto** for notebooks and reports  
- **Power BI** for dashboards  
- **GitHub Desktop** for version control  
- **VS Code** as the primary editor  

---

## Planned Workflow

1. Data ingestion & QC  
2. Shotgun preprocessing (milk microbiome)  
3. Alpha & beta diversity analysis  
4. Ecological type identification  
5. Longitudinal stability testing  
6. Visualization & dashboards  
7. Statistical modeling  
8. Integration with infant gut 16S data  
9. Interpretation & reporting  

---

## Ecological Types Identified (Milk Microbiome)

Using Bray–Curtis distances, PCoA, hierarchical clustering, and k-means silhouette scores, three ecological types were identified:

- **Ecotype 1 — Balanced community**  
  - Moderate richness  
  - Even distribution of common milk genera  
  - Represents the majority of participants  

- **Ecotype 2 — Streptococcus/Corynebacterium-dominant**  
  - Lower evenness  
  - One or two genera dominate  
  - Likely reflects maternal skin/oral transfer patterns  

- **Ecotype 3 — Outlier profile**  
  - Unique genus composition  
  - Driven by a single participant  
  - Important for understanding individual variability  

These ecological types are used in Script 03 to compare richness, evenness, and stability across lactation months.

---

## About This Project

This repository reflects my ongoing learning in microbiome analysis and early-life development. It is intentionally structured to mirror real-world workflows used in maternal/child health research, with clear separation between raw data, processed data, scripts, notebooks, and reports.

As the project evolves, each folder will be populated with reproducible code, visualizations, and documentation.

---

## Status

Early development: milk microbiome preprocessing, ecological typing, and diversity analysis completed.  
Next step: begin infant gut 16S MIMBES pipeline and integrate feeding-type metadata.


