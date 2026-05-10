rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)

load("phyloseqs/PS_merged_sanger.R")

MET <- as.data.frame(sample_data(PS_merged_sanger))

# --- Consolidate location names ---
# "Gagloaant ubani Tana" and "Kvelaantubani Tana" are the same locality as "Village Levitana"
MET$location[MET$location %in% c("Gagloaant ubani Tana", "Kvelaantubani Tana")] <- "Village Levitana"

cat("Location counts after renaming:\n")
print(table(MET$location))

# --- Create individual_ID from sample_ID ---
# sample_ID format: <individual>.<intestine_type> (e.g. 15.1 = upper, 15.2 = lower)
MET$individual_ID <- sub("\\..*", "", MET$sample_ID)

cat("\nindividual_ID counts:\n")
print(table(MET$individual_ID))

# --- Update phyloseq sample data ---
sample_data(PS_merged_sanger) <- sample_data(MET)

save(PS_merged_sanger, file = "phyloseqs/PS_merged_clean_sanger.R")
cat("\nPS_merged_clean_sanger saved.\n")
cat("Samples:", nsamples(PS_merged_sanger), "\n")
cat("ASVs:   ", ntaxa(PS_merged_sanger), "\n")
