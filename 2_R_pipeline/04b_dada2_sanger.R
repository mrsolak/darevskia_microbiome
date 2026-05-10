rm(list = ls())
setwd("D:/Rwd_claude/darevskia_16S_first")
library(dada2)

TRIM_DIR  <- "raw_files/2_trimmed_sanger"
FILT_DIR  <- "raw_files/3_filtered_sanger"
dir.create(FILT_DIR, showWarnings = FALSE)

# --- Input files ---
LIST    <- list.files(TRIM_DIR, full.names = TRUE)
F_reads <- LIST[grep("pair1.fastq.gz", LIST)]
R_reads <- LIST[grep("pair2.fastq.gz", LIST)]

sample.names <- gsub("-trimmed-pair1.fastq.gz", "", basename(F_reads))

filtFs <- file.path(FILT_DIR, paste0(sample.names, "_READ1_filt.fastq.gz"))
filtRs <- file.path(FILT_DIR, paste0(sample.names, "_READ2_filt.fastq.gz"))

# --- Quality filtering ---
track <- matrix(NA, nrow = length(sample.names), ncol = 2,
                dimnames = list(sample.names, c("input", "filtered")))

for (x in seq_along(F_reads)) {
  print(sample.names[x])
  out <- fastqPairedFilter(
    c(F_reads[x], R_reads[x]),
    c(filtFs[x], filtRs[x]),
    maxN = 0, maxEE = 2, minQ = 2, truncQ = 2,
    compress = TRUE, verbose = TRUE,
    minLen = c(270, 270), truncLen = c(270, 270)
  )
  track[x, "input"]    <- out["reads.in"]
  track[x, "filtered"] <- out["reads.out"]
}

# --- Dereplicate (only samples that passed filtering) ---
fnFs <- filtFs[file.exists(filtFs)]
fnRs <- filtRs[file.exists(filtRs)]
passed.names <- gsub("_READ1_filt.fastq.gz", "", basename(fnFs))

derepFs <- derepFastq(fnFs, n = 1e+05, verbose = TRUE)
derepRs <- derepFastq(fnRs, n = 1e+05, verbose = TRUE)
names(derepFs) <- passed.names
names(derepRs) <- passed.names

# --- Denoising ---
dadaFs <- dada(derepFs, err = NULL, multithread = FALSE, selfConsist = TRUE, MAX_CONSIST = 25)
dadaRs <- dada(derepRs, err = NULL, multithread = FALSE, selfConsist = TRUE, MAX_CONSIST = 25)

# --- Merging ---
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                      verbose = TRUE, minOverlap = 10, maxMismatch = 1,
                      justConcatenate = FALSE)

# --- Sequence table ---
seqtab <- makeSequenceTable(mergers)

# --- Read tracking ---
getN <- function(x) sum(getUniques(x))

track_full <- data.frame(
  sample     = sample.names,
  input      = track[sample.names, "input"],
  filtered   = track[sample.names, "filtered"],
  denoised_F = sapply(dadaFs[sample.names], getN),
  denoised_R = sapply(dadaRs[sample.names], getN),
  merged     = sapply(mergers[sample.names],  getN),
  row.names  = NULL
)

write.csv(track_full, "results/tables/dada2_read_tracking_sanger.csv", row.names = FALSE)
cat("Read tracking saved -> results/tables/dada2_read_tracking_sanger.csv\n")

# --- Save ---
save(seqtab, file = "phyloseqs/seqtab_sanger.R")
cat("seqtab saved -> phyloseqs/seqtab_sanger.R\n")
cat("Dimensions:", dim(seqtab), "\n")
