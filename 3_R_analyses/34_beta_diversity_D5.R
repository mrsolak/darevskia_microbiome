rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(vegan)
library(permute)
library(ape)
library(ggplot2)
library(dplyr)
library(patchwork)
library(flextable)
library(officer)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

# D5: D. portschinskii only — strata = individual_ID (paired design)
PS_D5 <- prune_samples(
  sample_data(PS_merged_sanger)$Species == "Darevskia portschinskii",
  PS_merged_sanger
)
PS_D5 <- prune_taxa(taxa_sums(PS_D5) > 0, PS_D5)
cat(sprintf("D5: %d samples, %d ASVs\n", nsamples(PS_D5), ntaxa(PS_D5)))

MET <- data.frame(as.data.frame(sample_data(PS_D5)), stringsAsFactors = FALSE)
MET$Life.stage <- factor(MET$Life.stage, levels = c("Juv", "Adult"))

dir.create("results/D5_portschinskii/plots/34_beta_diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D5_portschinskii/tables/34_beta_diversity", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1: Distance matrices
# ============================================================

cat("Computing distance matrices...\n")
dist_list <- list(
  "Bray-Curtis"        = phyloseq::distance(PS_D5, method = "bray"),
  "Jaccard"            = phyloseq::distance(PS_D5, method = "jaccard", binary = TRUE),
  "Weighted UniFrac"   = phyloseq::distance(PS_D5, method = "wunifrac"),
  "Unweighted UniFrac" = phyloseq::distance(PS_D5, method = "unifrac")
)
dist_names <- names(dist_list)

# ============================================================
# Step 2: PERMANOVA — within-individual strata (paired design)
# ============================================================

h_strata <- how(nperm = 999, blocks = MET$individual_ID)

cat("Running PERMANOVA (within-individual strata)...\n")
set.seed(42)
perm_strata <- lapply(dist_list, function(d)
  adonis2(d ~ Life.stage * sample_type, data = MET, permutations = h_strata, by = "terms"))

cat("\n=== PERMANOVA (within-individual strata) ===\n")
for (d in dist_names) { cat("\n---", d, "---\n"); print(perm_strata[[d]]) }

# ============================================================
# Step 3: PERMDISP — Life.stage and sample_type
# ============================================================

cat("\nRunning PERMDISP...\n")
set.seed(42)
permdisp_ls <- lapply(dist_list, function(d) {
  bd <- betadisper(d, MET$Life.stage)
  list(bd = bd, perm = permutest(bd, permutations = 999))
})
set.seed(42)
permdisp_type <- lapply(dist_list, function(d) {
  bd <- betadisper(d, MET$sample_type)
  list(bd = bd, perm = permutest(bd, permutations = 999))
})

cat("\n=== PERMDISP: Life stage ===\n")
for (d in dist_names) { cat("\n---", d, "---\n"); print(permdisp_ls[[d]]$perm) }

cat("\n=== PERMDISP: Sample type ===\n")
for (d in dist_names) { cat("\n---", d, "---\n"); print(permdisp_type[[d]]$perm) }

# ============================================================
# Step 4: PCoA plots
# ============================================================

# Life.stage colours; shape = sample_type
lifestage_colors <- c("Juv" = "#56B4E9", "Adult" = "#D55E00")
type_shapes      <- c("Upper intestine" = 16, "Lower intestine" = 17)

fmt_p_annot <- function(p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  p_str <- ifelse(p < 0.001, "<0.001", as.character(round(p, 3)))
  paste0(p_str, stars)
}

make_annot <- function(res) {
  # Terms may be Life.stage, sample_type, Life.stage:sample_type
  term_map <- list(
    "Life.stage"             = "Life stage",
    "sample_type"            = "Sample type",
    "Life.stage:sample_type" = "Interaction"
  )
  lines <- character(0)
  for (term in names(term_map)) {
    if (term %in% rownames(res)) {
      f  <- round(res[term, "F"],  2)
      r2 <- round(res[term, "R2"], 3)
      p  <- fmt_p_annot(res[term, "Pr(>F)"])
      lines <- c(lines, paste0(term_map[[term]], ": F=", f, " R\u00b2=", r2, " ", p))
    }
  }
  paste0("PERMANOVA (strata)\n", paste(lines, collapse = "\n"))
}

make_pcoa_plot <- function(dist_name, dist_obj, res) {
  pcoa_res <- ape::pcoa(dist_obj)
  scores   <- as.data.frame(pcoa_res$vectors[, 1:2])
  pct      <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)

  df <- data.frame(
    PC1        = scores[, 1],
    PC2        = scores[, 2],
    Life.stage = MET$Life.stage,
    sample_type = MET$sample_type,
    row.names  = rownames(scores)
  )

  annot <- make_annot(res)

  ggplot(df, aes(x = PC1, y = PC2, color = Life.stage)) +
    geom_point(aes(shape = sample_type), size = 2.8, alpha = 0.9) +
    stat_ellipse(aes(group = Life.stage), type = "t", level = 0.95, linewidth = 0.8) +
    scale_color_manual(values = lifestage_colors) +
    scale_shape_manual(values = type_shapes) +
    annotate("text", x = Inf, y = -Inf, label = annot,
             hjust = 1.05, vjust = -0.05, size = 2.7, lineheight = 1.3) +
    theme_bw(base_size = 12) +
    theme(legend.title  = element_text(face = "bold"),
          panel.grid    = element_blank(),
          plot.title    = element_text(face = "bold", size = 12)) +
    labs(title = dist_name,
         x     = paste0("PC1 (", pct[1], "%)"),
         y     = paste0("PC2 (", pct[2], "%)"),
         color = "Life stage",
         shape = "Sample type")
}

pcoa_plots <- mapply(
  make_pcoa_plot,
  dist_name = dist_names,
  dist_obj  = dist_list,
  res       = perm_strata,
  SIMPLIFY  = FALSE
)

p_combined <- wrap_plots(pcoa_plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("results/D5_portschinskii/plots/34_beta_diversity/PCoA_combined.png",
       p_combined, width = 13, height = 11, dpi = 300, bg = "white")

# ============================================================
# Step 5: Results tables as .docx
# ============================================================

fmt_p_tbl <- function(p) {
  if (is.na(p)) return("\u2014")
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  p_str <- ifelse(p < 0.001, "< 0.001", as.character(round(p, 3)))
  paste0(p_str, stars)
}

make_perm_table <- function(res, dist_name) {
  terms  <- c("Life.stage", "sample_type", "Life.stage:sample_type", "Residual", "Total")
  labels <- c("Life stage", "Sample type", "Life stage \u00d7 Sample type", "Residual", "Total")

  safe_val <- function(term, col) {
    if (term %in% rownames(res) && col %in% colnames(res)) res[term, col] else NA
  }

  df <- data.frame(
    Term = labels,
    df   = sapply(terms, function(t) safe_val(t, "Df")),
    F    = sapply(terms, function(t) { v <- safe_val(t, "F");  if (is.na(v)) "\u2014" else round(v, 3) }),
    R2   = sapply(terms, function(t) { v <- safe_val(t, "R2"); if (is.na(v)) "\u2014" else round(v, 3) }),
    p    = sapply(terms, function(t) fmt_p_tbl(safe_val(t, "Pr(>F)")))
  )
  colnames(df) <- c("Term", "df", "F", "R\u00b2", "p")

  flextable(df) %>%
    bold(part = "header") %>%
    set_caption(paste0("PERMANOVA (within-individual strata, sequential SS) \u2014 ", dist_name, " (D5)")) %>%
    footnote(i = 1, j = 5, part = "header",
             value = as_paragraph("Permutations restricted within individual_ID blocks (999 permutations)."),
             ref_symbols = "a") %>%
    autofit()
}

make_disp_table <- function(perm_ls, perm_type, dist_name) {
  extract_row <- function(res, label) {
    f <- round(res$perm$tab["Groups", "F"], 3)
    p <- res$perm$tab["Groups", "Pr(>F)"]
    data.frame(Group = label, F = f, p = fmt_p_tbl(p), stringsAsFactors = FALSE)
  }
  df <- rbind(extract_row(perm_ls,   "Life stage"),
              extract_row(perm_type, "Sample type"))
  colnames(df) <- c("Group", "F", "p")
  flextable(df) %>%
    bold(part = "header") %>%
    set_caption(paste0("PERMDISP \u2014 ", dist_name)) %>%
    autofit()
}

doc <- read_docx()
for (d in dist_names) {
  doc <- doc %>%
    body_add_par(d, style = "heading 1") %>%
    body_add_par("PERMANOVA", style = "heading 2") %>%
    body_add_flextable(make_perm_table(perm_strata[[d]], d)) %>%
    body_add_par("", style = "Normal") %>%
    body_add_par("PERMDISP (homogeneity of dispersions)", style = "heading 2") %>%
    body_add_flextable(make_disp_table(permdisp_ls[[d]], permdisp_type[[d]], d)) %>%
    body_add_par("", style = "Normal")
}
print(doc, target = "results/D5_portschinskii/tables/34_beta_diversity/beta_diversity_stats.docx")

cat("\nAll beta diversity outputs saved to:\n")
cat("  results/D5_portschinskii/plots/34_beta_diversity/\n")
cat("  results/D5_portschinskii/tables/34_beta_diversity/\n")
