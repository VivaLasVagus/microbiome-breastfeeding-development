###############################################
# 1. Load Metadata & File Paths
###############################################

meta <- read.csv("data/raw/infant_fastq/infant_metadata/metadata.csv", stringsAsFactors = FALSE)

fnFs <- setNames(meta$File1, meta$SampleID)
fnRs <- setNames(meta$File2, meta$SampleID)

meta$File1 <- file.path("data/raw/infant_fastq", meta$File1)
meta$File2 <- file.path("data/raw/infant_fastq", meta$File2)

fnFs <- setNames(meta$File1, meta$SampleID)
fnRs <- setNames(meta$File2, meta$SampleID)

head(fnFs)
head(fnRs)


###############################################
# 2. Load & Install DADA2
###############################################

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install dada2 if needed
if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2")
}

library(dada2)


###############################################
# 3. Create Filtered Output Directory
###############################################

if(!dir.exists("data/raw/infant_fastq/filtered")){
  dir.create("data/raw/infant_fastq/filtered")
}

filtFs <- file.path("data/raw/infant_fastq/filtered",
                    paste0(names(fnFs), "_F_filt.fastq.gz"))

filtRs <- file.path("data/raw/infant_fastq/filtered",
                    paste0(names(fnRs), "_R_filt.fastq.gz"))

packageVersion("dada2")


###############################################
# 4. Quality Profiles
###############################################

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


###############################################
# 5. Filtering & Trimming
###############################################

out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     truncLen=c(240,200),
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2,
                     rm.phix=TRUE,
                     compress=TRUE,
                     multithread=TRUE)

head(out)


###############################################
# 6. Dereplication
###############################################

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- names(filtFs)
names(derepRs) <- names(filtRs)


###############################################
# 7. Learn Error Rates
###############################################

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


###############################################
# 8. Sample Inference (DADA)
###############################################

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)


###############################################
# 9. Merge Paired Reads
###############################################

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

dadaFs[[1]]  # Inspect first sample


###############################################
# 10. Construct Sequence Table
###############################################

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


###############################################
# 11. Remove Chimeras
###############################################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)


###############################################
# 12. Assign Taxonomy (SILVA v138)
###############################################

taxa <- assignTaxonomy(
  seqtab.nochim,
  "data/silva/silva_nr_v138_train_set.fa.gz",
  multithread = FALSE
)

taxa <- addSpecies(
  taxa,
  "data/silva/silva_species_assignment_v138.fa.gz"
)
