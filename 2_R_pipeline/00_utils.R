# Custom functions by Jakub Kreisinger
# Ported from example_script_duplicates.Rmd

# Removes ASVs not present in both duplicates of each biological sample
dupl.concensus <- function(PHYLOSEQ, NAMES) {

  IDS      <- as.character(data.frame(sample_data(PHYLOSEQ))[, NAMES])
  IDS.dupl <- IDS[duplicated(IDS)]

  PHYLOSEQ <- prune_samples(IDS %in% IDS.dupl, PHYLOSEQ)
  if (length(IDS.dupl) * 2 < length(IDS)) {
    NONDUPLICATED <- prune_samples(!IDS %in% IDS.dupl, PHYLOSEQ)
    print(paste("Following names are nonduplicated", sample_names(NONDUPLICATED)))
  }

  CATS  <- as.character(data.frame(sample_data(PHYLOSEQ))[, NAMES])
  CATS2 <- levels(factor(CATS))
  OTU_TAB <- otu_table(PHYLOSEQ)
  rownames(OTU_TAB) <- CATS

  for (i in 1:length(CATS2)) {
    FILTER.act <- colSums(OTU_TAB[rownames(OTU_TAB) == CATS2[i], ] > 0) > 1
    OTU_TAB[rownames(OTU_TAB) == CATS2[i], ] <-
      t(apply(OTU_TAB[rownames(OTU_TAB) == CATS2[i], ], 1, function(x) x * FILTER.act))
  }

  rownames(OTU_TAB) <- sample_names(PHYLOSEQ)
  otu_table(PHYLOSEQ) <- OTU_TAB
  prune_taxa(taxa_sums(PHYLOSEQ) > 0, PHYLOSEQ)
}

# Merges duplicate sequencing samples into one sample per biological ID
merge.duplicates <- function(PHYLOSEQ, NAMES) {
  CATS <- as.character(data.frame(sample_data(PHYLOSEQ))[, NAMES])
  sample_data(PHYLOSEQ)$duplic.id <- CATS
  SAMDAT     <- sample_data(PHYLOSEQ)
  SAMDAT.sub <- SAMDAT[duplicated(CATS) == FALSE, ]
  FASTA      <- refseq(PHYLOSEQ)
  rownames(SAMDAT.sub) <- SAMDAT.sub$duplic.id
  PHYLOSEQ.merge <- merge_samples(PHYLOSEQ, "duplic.id")
  sample_data(PHYLOSEQ.merge) <- SAMDAT.sub
  merge_phyloseq(PHYLOSEQ.merge, FASTA)
}
