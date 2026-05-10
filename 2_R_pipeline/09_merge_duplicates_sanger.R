rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(ggplot2)
library(vegan)

source("scripts/2_R_scripts/00_utils.R")

load("phyloseqs/PS_meta_sanger.R")

dir.create("results/tables/removed_reports",    recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots/04_duplicates_sanger", recursive = TRUE, showWarnings = FALSE)

MET     <- as.data.frame(sample_data(PS_meta_sanger))
SUMMARY <- summary(as.factor(MET$sample_ID), maxsum = 10000)

cat("Duplicate counts per biological sample:\n")
print(rev(sort(SUMMARY)))

# --- Check for triplicates / tetraplicates ---
high_rep <- SUMMARY[SUMMARY >= 3]
if (length(high_rep) > 0) {
  write.csv(as.data.frame(high_rep),
            "results/tables/removed_reports/triplicates_tetraplicates_sanger.csv",
            row.names = FALSE)
  cat("WARNING: triplicates or tetraplicates found:\n")
  print(high_rep)
  stop("Resolve triplicates/tetraplicates before continuing.")
}
cat("No triplicates or tetraplicates found.\n")

# --- Report and remove singletons ---
singleton_ids   <- names(SUMMARY[SUMMARY == 1])
singleton_samps <- rownames(MET)[MET$sample_ID %in% singleton_ids]
cat("\nSingletons found:", length(singleton_samps), "\n")

if (length(singleton_samps) > 0) {
  print(singleton_samps)
  write.csv(
    data.frame(sample    = singleton_samps,
               sample_ID = MET[singleton_samps, "sample_ID"],
               reads     = sample_sums(PS_meta_sanger)[singleton_samps]),
    "results/tables/removed_reports/removed_singletons_sanger.csv",
    row.names = FALSE
  )
} else {
  write.csv(data.frame(sample = character(0), sample_ID = character(0), reads = numeric(0)),
            "results/tables/removed_reports/removed_singletons_sanger.csv",
            row.names = FALSE)
  cat("No singletons — all samples are complete pairs.\n")
}

SUMMARY.n2     <- names(SUMMARY[SUMMARY == 2])
PS_dupl_sanger <- prune_samples(MET$sample_ID %in% SUMMARY.n2, PS_meta_sanger)
PS_dupl_sanger <- prune_taxa(taxa_sums(PS_dupl_sanger) > 0, PS_dupl_sanger)

cat("\nSamples remaining (complete pairs):", nsamples(PS_dupl_sanger), "\n")
cat("Biological samples:", nsamples(PS_dupl_sanger) / 2, "\n")

save(PS_dupl_sanger, file = "phyloseqs/PS_dupl_sanger.R")

# --- Procrustes: compare duplicate sets ---
MET_d  <- as.data.frame(sample_data(PS_dupl_sanger))
dupl.1 <- prune_samples(duplicated(MET_d$sample_ID) == FALSE, PS_dupl_sanger)
dupl.2 <- prune_samples(duplicated(MET_d$sample_ID) == TRUE,  PS_dupl_sanger)

BC1 <- as.matrix(vegdist(otu_table(transform_sample_counts(dupl.1, function(x) x / sum(x)))))
BC2 <- as.matrix(vegdist(otu_table(transform_sample_counts(dupl.2, function(x) x / sum(x)))))

rownames(BC1) <- colnames(BC1) <- sample_data(dupl.1)$sample_ID
rownames(BC2) <- colnames(BC2) <- sample_data(dupl.2)$sample_ID
BC2 <- BC2[rownames(BC1), rownames(BC1)]

BC1.pcoa <- cmdscale(as.dist(BC1))
BC2.pcoa <- cmdscale(as.dist(BC2))

PROTEST <- protest(BC1.pcoa, BC2.pcoa)
cat("\n=== Protest result ===\n")
print(PROTEST)

write.csv(
  data.frame(correlation  = PROTEST$t0,
             p_value      = PROTEST$signif,
             permutations = 999),
  "results/tables/removed_reports/protest_result_sanger.csv",
  row.names = FALSE
)

DF.prot <- data.frame(PROTEST$Y, PROTEST$X, sample_data(dupl.1),
                      resid = resid(PROTEST))
names(DF.prot)[1:4] <- c("X1", "X2", "Y1", "Y2")

p_protest <- ggplot(DF.prot, aes(x = X1, y = X2)) +
  geom_point(size = 2) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2),
               arrow = arrow(length = unit(0.1, "cm"))) +
  theme_bw(base_size = 12) +
  labs(title    = paste0("Procrustes: correlation = ", round(PROTEST$t0, 3),
                         ", p = ", PROTEST$signif),
       subtitle = "Arrow length = distance between duplicates",
       x = "PCoA axis 1", y = "PCoA axis 2")

ggsave("results/plots/04_duplicates_sanger/procrustes_duplicates_sanger.png",
       p_protest, width = 7, height = 6, dpi = 300, bg = "white")
cat("Procrustes plot saved.\n")

# --- Consensus filter + merge ---
snapshot <- function(ps) {
  c(samples = nsamples(ps), asvs = ntaxa(ps), reads = sum(sample_sums(ps)))
}

s_before  <- snapshot(PS_dupl_sanger)
PS_consist <- dupl.concensus(PHYLOSEQ = PS_dupl_sanger, NAMES = "sample_ID")
s_consist  <- snapshot(PS_consist)

PS_merged_sanger <- merge.duplicates(PHYLOSEQ = PS_consist, NAMES = "sample_ID")
s_merge          <- snapshot(PS_merged_sanger)

# --- Remove samples below 500 reads post-merge ---
low_depth <- sample_sums(PS_merged_sanger)[sample_sums(PS_merged_sanger) < 500]
write.csv(
  data.frame(sample = names(low_depth), reads = low_depth),
  "results/tables/removed_reports/removed_low_depth_post_merge_sanger.csv",
  row.names = FALSE
)
cat("\nSamples below 500 reads after merging:", length(low_depth), "\n")
if (length(low_depth) > 0) print(sort(low_depth))

PS_merged_sanger <- prune_samples(sample_sums(PS_merged_sanger) >= 500, PS_merged_sanger)
PS_merged_sanger <- prune_taxa(taxa_sums(PS_merged_sanger) > 0, PS_merged_sanger)
s_final          <- snapshot(PS_merged_sanger)

# --- Before/after summary ---
report <- data.frame(
  step    = c("PS_dupl_sanger (input)",
              "after dupl.concensus",
              "after merge.duplicates",
              "after low-depth filter + zero ASV drop"),
  samples = c(s_before["samples"], s_consist["samples"],
              s_merge["samples"],  s_final["samples"]),
  asvs    = c(s_before["asvs"],    s_consist["asvs"],
              s_merge["asvs"],     s_final["asvs"]),
  reads   = c(s_before["reads"],   s_consist["reads"],
              s_merge["reads"],    s_final["reads"])
)
write.csv(report,
          "results/tables/removed_reports/merge_duplicates_summary_sanger.csv",
          row.names = FALSE)
cat("\n=== Merge duplicates summary ===\n")
print(report)

save(PS_merged_sanger, file = "phyloseqs/PS_merged_sanger.R")
cat("\nPS_merged_sanger saved.\n")
