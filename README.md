# microbiome-breastfeeding-development

A reproducible analysis exploring how breastfeeding shapes early-life microbiome ecology — beginning with the human milk microbiome — and how these microbial patterns may relate to infant gut development. This project is part of my learning and portfolio-building journey in microbiome science, maternal/child health, and early-life development.

---

## Project Goals

- Characterize the human milk microbiome using shotgun-derived genus-level data  
- Identify ecological types (community structure patterns) in milk samples  
- Examine richness, evenness, and longitudinal stability across lactation months  
- Build reproducible R workflows for microbiome analysis (phyloseq, vegan, tidyverse)  
- Integrate infant gut 16S rRNA data using a reproducible DADA2 pipeline
- Compare milk microbiome ecological types with infant gut community development
- Build longitudinal models linking microbial features to feeding type and developmental outcomes 
- Practice reproducible workflows using GitHub, R, Quarto, and tidy data principles
- Practice work in QIIME2

---

## Repository Structure

data/
raw/            # Unmodified source data
dryad_milk_shotgun/          # Dryad data
infant_fastq/                # INFANTMET cohort data, Raw FASTQ files from SRR accessions
infant_metadata/             # metadata.csv linking Infant ID, Week, SRR, FASTQ filenames
processed/                   # Cleaned and transformed phyloseq objects

scripts/
01_qc_preprocessing.R        # Data ingestion, filtering, genus-level collapse
02_ecological_typing.R       # Bray–Curtis PCoA, clustering, silhouette scores
03_ecological_diversity.R    # Richness, evenness, stability, ecotype comparisons

**Currently Working**: 
04_dada2_pipeline.R          # Infant gut 16S DADA2 workflow
05_qiime2_aside.R            # QIIME2 learning module and comparison

notebooks/
exploratory_analysis.qmd     # Quarto notebooks for exploratory work
functional_pathways.qmd

reports/
microbiome_summary.md        # Narrative summaries
developmental_outcomes.md
summary/
ecological_types.md          # Ecotype definitions & results summary

dashboards/
powerbi/                     # Power BI files for interactive visualization

results/                     #Plots and .csv files from R


---

## Tools & Technologies

- **R** (tidyverse, phyloseq, vegan, microbiome, ggplot2)  
- **Quarto** for notebooks and reports  
- **Power BI** for dashboards  
- **GitHub Desktop** for version control  
- **RStudio** as primary editor

---

## Planned Workflow

**Milk Microbiome Workflow**
1. Data ingestion & QC  
2. Shotgun preprocessing (milk microbiome)  
3. Alpha & beta diversity analysis  
4. Ecological type identification  
5. Longitudinal stability testing  
6. Visualization & dashboards  
7. Statistical modeling  
8. Integration with infant gut 16S data  
9. Interpretation & reporting

**Infant Gut Workflow**
10. Download SRR accessions and convert .sra → paired FASTQs
11. Build metadata table linking InfantID, Week, SRR, and FASTQ files
12. Run full DADA2 pipeline: filtering, dereplication, error learning, denoising, merging, chimera removal
13. Construct ASV table and assign taxonomy
14. Build phyloseq object for infant gut samples
15. Compare infant gut community development across weeks
16. Optional: run parallel QIIME 2 workflow for learning and cross-platform comparison

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
## Infant Gut 16S rRNA Pipeline

Project includes a full infant gut microbiome workflow using paired-end 16S rRNA sequencing data from 10 infants sampled at Week 1, Week 8, and Week 24.

### Cohort Description (INFANTMET)
- 10 infants
- 3 timepoints: Week 1, Week 8, Week 24
- 30 paired-end 16S rRNA samples
- SRR accessions sourced from NCBI SRA
- FASTQs generated using fastq-dump --split-files --gzip


Data Processing Steps
- Download .sra files via NCBI SRA Toolkit
- Convert .sra → paired FASTQs using fastq-dump --split-files --gzip
- Build metadata table linking InfantID, Week, SRR, and FASTQ filenames
- Run full DADA2 pipeline (v1.40.0), including:
    - Filtering & trimming
    - Dereplication
    - Error model learning
    - Denoising
    - Merging paired reads
    - Chimera removal
    - ASV table construction
    - Taxonomy assignment
- Build phyloseq object for downstream analysis
- Compare infant gut trajectories across early life

Optional QIIME 2 Aside
- A parallel QIIME 2 workflow is included for learning and cross-platform comparison.

---

## About This Project

This repository reflects my ongoing learning in microbiome analysis and early-life development. It is intentionally structured to mirror real-world workflows used in maternal/child health research, with clear separation between raw data, processed data, scripts, notebooks, and reports.

As the project evolves, each folder will be populated with reproducible code, visualizations, and documentation.

---

## Status

Milk microbiome preprocessing, ecological typing, and diversity analysis completed.
Infant gut 16S pipeline initialized:

✔ SRR dataset identified

✔ .sra files downloaded

✔ FASTQs generated

✔ Metadata table built

✔ DADA2 installed and configured

✔ Filtering step ready to run

✔ Metadata successfully linked to FASTQ file paths in R

Next step: run DADA2 filtering and begin infant gut ASV construction.

