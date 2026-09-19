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
# Create filtered outputs
###############################
BiocManager::install("dada2")

library(dada2)

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


