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

# D3: upper intestine only
PS_D3 <- prune_samples(
  sample_data(PS_merged_sanger)$sample_type == "Upper intestine",
  PS_merged_sanger
)
PS_D3 <- prune_taxa(taxa_sums(PS_D3) > 0, PS_D3)
cat(sprintf("D3: %d samples, %d ASVs\n", nsamples(PS_D3), ntaxa(PS_D3)))

MET <- data.frame(as.data.frame(sample_data(PS_D3)), stringsAsFactors = FALSE)
OTU <- as.data.frame(as.matrix(otu_table(PS_D3)))
if (!taxa_are_rows(PS_D3)) OTU <- t(OTU)

dir.create("results/D3_upper/plots/25_venn_diagram", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D3_upper/tables/25_venn_diagram", recursive = TRUE, showWarnings = FALSE)

get_asv_set <- function(species, prev_threshold = 0) {
  samps <- rownames(MET)[MET$Species == species]
  n     <- length(samps)
  prev  <- rowSums(OTU[, samps, drop = FALSE] > 0)
  min_samples <- max(1L, ceiling(prev_threshold * n))
  names(prev)[prev >= min_samples]
}

# 2 groups: DD Upper vs DP Upper
group_labels  <- c("DD Upper", "DP Upper")
group_species <- c("Darevskia dahli", "Darevskia portschinskii")

thresholds <- c(0, 0.25, 0.50)
threshold_labels <- c(
  "No threshold\n(\u2265 1 sample)",
  "\u2265 25% of samples",
  "\u2265 50% of samples\n(core microbiome)"
)

cat("=== Sample counts per group (D3) ===\n")
for (i in seq_along(group_labels)) {
  samps <- rownames(MET)[MET$Species == group_species[i]]
  cat(sprintf("  %-22s  n = %d\n", group_labels[i], length(samps)))
}

all_sets <- lapply(thresholds, function(thr) {
  setNames(
    lapply(group_species, function(sp) get_asv_set(sp, thr)),
    group_labels
  )
})

for (k in seq_along(thresholds)) {
  cat(sprintf("\n=== Threshold: %s ===\n", threshold_labels[k]))
  for (nm in group_labels) {
    cat(sprintf("  %-22s  %d ASVs\n", nm, length(all_sets[[k]][[nm]])))
  }
  cat("  Total unique:", length(unique(unlist(all_sets[[k]]))), "\n")
}

make_venn2 <- function(sets, title) {
  p <- ggVennDiagram(
    sets,
    label       = "count",
    label_alpha = 0,
    edge_size   = 0.8
  )

  # Push set labels outward from centroid
  if (length(p$layers) >= 3) {
    lbl <- p$layers[[3]]$data
    cx  <- mean(lbl$X)
    cy  <- mean(lbl$Y)
    p$layers[[3]]$data$X <- cx + (lbl$X - cx) * 1.5
    p$layers[[3]]$data$Y <- cy + (lbl$Y - cy) * 1.5
    p$coordinates$clip <- "off"
  }

  p +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      name      = "ASVs\nin region"
    ) +
    scale_color_manual(values = c("#E69F00", "#0072B2")) +
    labs(title = title) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(25, 35, 25, 35)
    )
}

panels <- mapply(make_venn2, sets = all_sets, title = threshold_labels, SIMPLIFY = FALSE)

p_combined <- wrap_plots(panels, nrow = 1) +
  plot_annotation(
    title = "ASV overlap: D. dahli vs D. portschinskii — D3 upper intestine",
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
  )

ggsave("results/D3_upper/plots/25_venn_diagram/venn_2way_thresholds.png",
       p_combined, width = 14, height = 6, dpi = 300, bg = "white")

p_single <- make_venn2(all_sets[[1]], "No threshold (\u2265 1 sample)")
ggsave("results/D3_upper/plots/25_venn_diagram/venn_2way.png",
       p_single, width = 7, height = 6, dpi = 300, bg = "white")

set_names_short <- c("dahli_upper", "port_upper")

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
            sprintf("results/D3_upper/tables/25_venn_diagram/venn_summary_thr%s.csv", thr_tag),
            row.names = FALSE)
  write.csv(as.data.frame(overlap_mat),
            sprintf("results/D3_upper/tables/25_venn_diagram/venn_overlaps_thr%s.csv", thr_tag),
            row.names = TRUE)
}

cat("\nOutputs saved to:\n")
cat("  results/D3_upper/plots/25_venn_diagram/\n")
cat("  results/D3_upper/tables/25_venn_diagram/\n")
