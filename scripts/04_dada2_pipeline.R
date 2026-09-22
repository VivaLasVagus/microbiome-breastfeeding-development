meta <- read.csv("data/raw/infant_fastq/infant_metadata/metadata.csv", stringsAsFactors = FALSE)

fnFs <- setNames(meta$File1, meta$SampleID)
fnRs <- setNames(meta$File2, meta$SampleID)


meta$File1 <- file.path("data/raw/infant_fastq", meta$File1)
meta$File2 <- file.path("data/raw/infant_fastq", meta$File2)

fnFs <- setNames(meta$File1, meta$SampleID)
fnRs <- setNames(meta$File2, meta$SampleID)

head(fnFs)
head(fnRs)

###############################
# Set up DADA2 package
###############################

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install dada2 if needed
if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2")
}

# Load dada2
library(dada2)


###############################
# Create filtered outputs
###############################

# Create filtered directory if it doesn't exist
if(!dir.exists("data/raw/infant_fastq/filtered")){
  dir.create("data/raw/infant_fastq/filtered")
}

filtFs <- file.path("data/raw/infant_fastq/filtered",
                    paste0(names(fnFs), "_F_filt.fastq.gz"))

filtRs <- file.path("data/raw/infant_fastq/filtered",
                    paste0(names(fnRs), "_R_filt.fastq.gz"))

packageVersion("dada2")


plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

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

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- names(filtFs)
names(derepRs) <- names(filtRs)

