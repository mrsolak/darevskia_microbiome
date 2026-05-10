rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)

load("phyloseqs/PS_raw_sanger.R")

MET    <- read.csv("phyloseqs/Darevskia_16S_merged_mari.csv")
SNAMES <- sample_names(PS_raw_sanger)
SN     <- MET[["Primers_Plate_Indexes_Seqs_miseq_orient"]]

dir.create("results/tables/removed_reports", recursive = TRUE, showWarnings = FALSE)

# --- PS samples missing from metadata (will be removed) ---
ps_not_in_met <- setdiff(SNAMES, SN)
write.csv(
  data.frame(sample = ps_not_in_met,
             reads  = sample_sums(PS_raw_sanger)[ps_not_in_met]),
  "results/tables/removed_reports/removed_PS_not_in_metadata_sanger.csv",
  row.names = FALSE
)
cat("PS samples not in metadata (removed):", length(ps_not_in_met), "\n")
print(ps_not_in_met)

# --- Metadata rows with no PS sample (reported only, not removed) ---
met_not_in_ps <- setdiff(SN, SNAMES)
write.csv(
  MET[MET[["Primers_Plate_Indexes_Seqs_miseq_orient"]] %in% met_not_in_ps, ],
  "results/tables/removed_reports/metadata_rows_not_in_PS_sanger.csv",
  row.names = FALSE
)
cat("Metadata rows with no sequencing data (reported only):", length(met_not_in_ps), "\n")
print(met_not_in_ps)

# --- Remove PS samples missing from metadata ---
PS_raw_sanger <- prune_samples(!(SNAMES %in% ps_not_in_met), PS_raw_sanger)

# --- Merge metadata ---
rownames(MET) <- MET[["Primers_Plate_Indexes_Seqs_miseq_orient"]]
PS_meta_sanger <- merge_phyloseq(PS_raw_sanger, sample_data(MET))

# --- Report and remove controls ---
is_control <- grepl("^NK_|^PK_", sample_data(PS_meta_sanger)[["sample_ID"]])
controls_df <- data.frame(
  sample    = sample_names(PS_meta_sanger)[is_control],
  sample_ID = sample_data(PS_meta_sanger)[["sample_ID"]][is_control],
  reads     = sample_sums(PS_meta_sanger)[is_control]
)
write.csv(controls_df,
          "results/tables/removed_reports/removed_controls_sanger.csv",
          row.names = FALSE)
cat("\nControls found and removed:", sum(is_control), "\n")
print(controls_df)

PS_meta_sanger <- prune_samples(!is_control, PS_meta_sanger)
PS_meta_sanger <- prune_taxa(taxa_sums(PS_meta_sanger) > 0, PS_meta_sanger)

save(PS_meta_sanger, file = "phyloseqs/PS_meta_sanger.R")

cat("\n=== PS_meta_sanger summary (controls removed) ===\n")
cat("Samples:", nsamples(PS_meta_sanger), "\n")
cat("ASVs:   ", ntaxa(PS_meta_sanger), "\n")
