rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(ShortRead)

load("phyloseqs/seqtab_sanger.R")    # loads: seqtab
load("phyloseqs/tax_table_sanger.R") # loads: taxa

OTU <- otu_table(seqtab, taxa_are_rows = FALSE)
TAX <- tax_table(taxa)
ASV <- readDNAStringSet("phyloseqs/non_chimeras_sanger.fasta")
names(ASV) <- as.character(ASV)  # OTU/TAX tables use sequences as names, not ASV IDs

PS_raw_sanger <- merge_phyloseq(OTU, TAX, ASV)

dir.create("results/tables/removed_reports", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/sample_sums",     recursive = TRUE, showWarnings = FALSE)

# --- Remove ASVs with NA at Phylum level ---
na_phylum <- is.na(as.vector(tax_table(PS_raw_sanger)[, "Phylum"]))
removed_na <- data.frame(
  ASV       = taxa_names(PS_raw_sanger)[na_phylum],
  abundance = taxa_sums(PS_raw_sanger)[na_phylum],
  as.data.frame(tax_table(PS_raw_sanger)[na_phylum, ])
)
write.csv(removed_na,
          "results/tables/removed_reports/removed_NA_phylum_sanger.csv",
          row.names = FALSE)
cat("Removed (NA Phylum):", sum(na_phylum), "ASVs,",
    sum(removed_na$abundance), "reads\n")

PS_raw_sanger <- prune_taxa(!na_phylum, PS_raw_sanger)

# --- Remove Chloroplast-order ASVs ---
orders         <- as.vector(tax_table(PS_raw_sanger)[, "Order"])
is_chloroplast <- !is.na(orders) & orders == "Chloroplast"
removed_chloro <- data.frame(
  ASV       = taxa_names(PS_raw_sanger)[is_chloroplast],
  abundance = taxa_sums(PS_raw_sanger)[is_chloroplast],
  as.data.frame(tax_table(PS_raw_sanger)[is_chloroplast, ])
)
write.csv(removed_chloro,
          "results/tables/removed_reports/removed_Chloroplast_sanger.csv",
          row.names = FALSE)
cat("Removed (Chloroplast):", sum(is_chloroplast), "ASVs,",
    sum(removed_chloro$abundance), "reads\n")

PS_raw_sanger <- prune_taxa(!is_chloroplast, PS_raw_sanger)

# --- Save ---
save(PS_raw_sanger, file = "phyloseqs/PS_raw_sanger.R")

# --- Report sample read depths ---
sums_df <- data.frame(
  sample = sample_names(PS_raw_sanger),
  reads  = sample_sums(PS_raw_sanger)
)
write.csv(sums_df,
          "results/tables/sample_sums/sample_sums_PS_raw_sanger.csv",
          row.names = FALSE)

cat("\n=== PS_raw_sanger summary ===\n")
cat("Samples:", nsamples(PS_raw_sanger), "\n")
cat("ASVs:   ", ntaxa(PS_raw_sanger), "\n")
cat("Reads:  ", sum(sums_df$reads), "\n\n")
print(summary(sums_df$reads))
