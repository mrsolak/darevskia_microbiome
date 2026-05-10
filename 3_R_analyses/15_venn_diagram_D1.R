rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(patchwork)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

dir.create("results/D1_all/plots/15_venn_diagram", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D1_all/tables/15_venn_diagram", recursive = TRUE, showWarnings = FALSE)

MET <- data.frame(as.data.frame(sample_data(PS_merged_sanger)), stringsAsFactors = FALSE)
OTU <- as.data.frame(as.matrix(otu_table(PS_merged_sanger)))
if (!taxa_are_rows(PS_merged_sanger)) OTU <- t(OTU)
# OTU is now ASVs × samples

# ============================================================
# Helper: ASV sets at a given prevalence threshold
# prev_threshold = minimum proportion of samples in the group
# that must have > 0 reads for an ASV to be included
# ============================================================

get_asv_set <- function(species, type, prev_threshold = 0) {
  samps <- rownames(MET)[MET$Species == species & MET$sample_type == type]
  n     <- length(samps)
  # count how many samples have > 0 reads for each ASV
  prev  <- rowSums(OTU[, samps, drop = FALSE] > 0)
  min_samples <- max(1L, ceiling(prev_threshold * n))
  names(prev)[prev >= min_samples]
}

group_labels <- c("DD Upper", "DD Lower", "DP Upper", "DP Lower")
group_species <- c("Darevskia dahli", "Darevskia dahli",
                   "Darevskia portschinskii", "Darevskia portschinskii")
group_types   <- c("Upper intestine", "Lower intestine",
                   "Upper intestine", "Lower intestine")

thresholds <- c(0, 0.25, 0.50)
threshold_labels <- c(
  "No threshold\n(≥ 1 sample)",
  "≥ 25% of samples",
  "≥ 50% of samples\n(core microbiome)"
)

# ============================================================
# Print per-group sample counts
# ============================================================

cat("=== Sample counts per group ===\n")
for (i in seq_along(group_labels)) {
  samps <- rownames(MET)[MET$Species == group_species[i] & MET$sample_type == group_types[i]]
  cat(sprintf("  %-22s  n = %d\n", gsub("\n", " ", group_labels[i]), length(samps)))
}

# ============================================================
# Build sets and print summaries for each threshold
# ============================================================

all_sets <- lapply(thresholds, function(thr) {
  s <- setNames(
    mapply(get_asv_set, group_species, group_types,
           MoreArgs = list(prev_threshold = thr), SIMPLIFY = FALSE),
    group_labels
  )
  s
})

for (k in seq_along(thresholds)) {
  cat(sprintf("\n=== Threshold: %s ===\n", threshold_labels[k]))
  for (nm in group_labels) {
    cat(sprintf("  %-22s  %d ASVs\n", gsub("\n", " ", nm), length(all_sets[[k]][[nm]])))
  }
  cat("  Total unique:", length(unique(unlist(all_sets[[k]]))), "\n")
}

# ============================================================
# Venn diagram function
# ============================================================

make_venn <- function(sets, title) {
  p <- ggVennDiagram(
    sets,
    label       = "count",
    label_alpha = 0,
    edge_size   = 0.8
  )

  # Layer 3 (GeomText, 4 rows) holds the set labels — push them away from centre
  lbl <- p$layers[[3]]$data
  cx  <- mean(lbl$X)
  cy  <- mean(lbl$Y)
  p$layers[[3]]$data$X <- cx + (lbl$X - cx) * 1.55
  p$layers[[3]]$data$Y <- cy + (lbl$Y - cy) * 1.55

  # Disable clipping so pushed labels are not cut off at panel edge
  p$coordinates$clip <- "off"

  p +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      name      = "ASVs\nin region"
    ) +
    scale_color_manual(values = c("#E69F00", "#E69F00", "#0072B2", "#0072B2")) +
    labs(title = title) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(25, 35, 25, 35)
    )
}

panels <- mapply(make_venn, sets = all_sets, title = threshold_labels, SIMPLIFY = FALSE)

p_combined <- wrap_plots(panels, nrow = 1) +
  plot_annotation(
    title = "ASV overlap across species and intestine types",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )
  )

ggsave("results/D1_all/plots/15_venn_diagram/venn_4way_thresholds.png",
       p_combined, width = 18, height = 7, dpi = 300, bg = "white")

# keep the original single-panel file too
p_single <- make_venn(all_sets[[1]], "No threshold (≥ 1 sample)")
ggsave("results/D1_all/plots/15_venn_diagram/venn_4way.png",
       p_single, width = 8, height = 7, dpi = 300, bg = "white")

# ============================================================
# Summary tables — one per threshold
# ============================================================

set_names_short <- c("dahli_upper", "dahli_lower", "port_upper", "port_lower")

for (k in seq_along(thresholds)) {
  sets <- all_sets[[k]]
  thr_tag <- gsub("[^0-9]", "", as.character(thresholds[k]))
  if (thr_tag == "") thr_tag <- "0"

  n <- length(sets)
  overlap_mat <- matrix(NA_integer_, n, n,
                        dimnames = list(set_names_short, set_names_short))
  for (i in seq_len(n))
    for (j in seq_len(n))
      overlap_mat[i, j] <- length(intersect(sets[[i]], sets[[j]]))

  unique_counts <- sapply(seq_len(n), function(i) {
    others <- unlist(sets[-i])
    sum(!sets[[i]] %in% others)
  })

  summary_df <- data.frame(
    Group       = group_labels,
    Total_ASVs  = diag(overlap_mat),
    Unique_ASVs = unique_counts,
    stringsAsFactors = FALSE
  )

  write.csv(summary_df,
            sprintf("results/D1_all/tables/15_venn_diagram/venn_summary_thr%s.csv", thr_tag),
            row.names = FALSE)
  write.csv(as.data.frame(overlap_mat),
            sprintf("results/D1_all/tables/15_venn_diagram/venn_overlaps_thr%s.csv", thr_tag),
            row.names = TRUE)
}

cat("\nOutputs saved to:\n")
cat("  results/D1_all/plots/15_venn_diagram/venn_4way_thresholds.png  (3-panel combined)\n")
cat("  results/D1_all/plots/15_venn_diagram/venn_4way.png             (no-threshold single panel)\n")
cat("  results/D1_all/tables/15_venn_diagram/\n")
