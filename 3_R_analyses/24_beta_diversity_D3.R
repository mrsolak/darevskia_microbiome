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

# D3: upper intestine only — free permutation (each individual appears once)
PS_D3 <- prune_samples(
  sample_data(PS_merged_sanger)$sample_type == "Upper intestine",
  PS_merged_sanger
)
PS_D3 <- prune_taxa(taxa_sums(PS_D3) > 0, PS_D3)
cat(sprintf("D3: %d samples, %d ASVs\n", nsamples(PS_D3), ntaxa(PS_D3)))

MET <- data.frame(as.data.frame(sample_data(PS_D3)), stringsAsFactors = FALSE)

dir.create("results/D3_upper/plots/24_beta_diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D3_upper/tables/24_beta_diversity", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1: Distance matrices
# ============================================================

cat("Computing distance matrices...\n")
dist_list <- list(
  "Bray-Curtis"        = phyloseq::distance(PS_D3, method = "bray"),
  "Jaccard"            = phyloseq::distance(PS_D3, method = "jaccard", binary = TRUE),
  "Weighted UniFrac"   = phyloseq::distance(PS_D3, method = "wunifrac"),
  "Unweighted UniFrac" = phyloseq::distance(PS_D3, method = "unifrac")
)
dist_names <- names(dist_list)

# ============================================================
# Step 2: PERMANOVA — free permutation (unpaired, each individual once)
# ============================================================

h_free <- how(nperm = 999)

cat("Running PERMANOVA (free permutation)...\n")
set.seed(42)
perm_free <- lapply(dist_list, function(d)
  adonis2(d ~ Species, data = MET, permutations = h_free, by = "terms"))

cat("\n=== PERMANOVA (free permutation) ===\n")
for (d in dist_names) { cat("\n---", d, "---\n"); print(perm_free[[d]]) }

# ============================================================
# Step 3: PERMDISP — Species only
# ============================================================

cat("\nRunning PERMDISP...\n")
set.seed(42)
permdisp_sp <- lapply(dist_list, function(d) {
  bd <- betadisper(d, MET$Species)
  list(bd = bd, perm = permutest(bd, permutations = 999))
})

cat("\n=== PERMDISP: Species ===\n")
for (d in dist_names) { cat("\n---", d, "---\n"); print(permdisp_sp[[d]]$perm) }

# ============================================================
# Step 4: PCoA plots (Species only)
# ============================================================

species_colors <- c("Darevskia dahli" = "#E69F00", "Darevskia portschinskii" = "#0072B2")

fmt_p_annot <- function(p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  p_str <- ifelse(p < 0.001, "<0.001", as.character(round(p, 3)))
  paste0(p_str, stars)
}

make_annot <- function(res) {
  f  <- round(res["Species", "F"],  2)
  r2 <- round(res["Species", "R2"], 3)
  p  <- fmt_p_annot(res["Species", "Pr(>F)"])
  paste0("PERMANOVA (free)\nSpecies: F=", f, " R\u00b2=", r2, " ", p)
}

make_pcoa_plot <- function(dist_name, dist_obj, res) {
  pcoa_res <- ape::pcoa(dist_obj)
  scores   <- as.data.frame(pcoa_res$vectors[, 1:2])
  pct      <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)

  df <- data.frame(
    PC1     = scores[, 1],
    PC2     = scores[, 2],
    Species = MET$Species,
    row.names = rownames(scores)
  )

  annot <- make_annot(res)

  ggplot(df, aes(x = PC1, y = PC2, color = Species)) +
    geom_point(size = 2.8, alpha = 0.9) +
    stat_ellipse(aes(group = Species), type = "t", level = 0.95, linewidth = 0.8) +
    scale_color_manual(values = species_colors,
                       labels = c("Darevskia dahli"         = "D. dahli",
                                  "Darevskia portschinskii" = "D. portschinskii")) +
    annotate("text", x = Inf, y = -Inf, label = annot,
             hjust = 1.05, vjust = -0.05, size = 2.9, lineheight = 1.3) +
    theme_bw(base_size = 12) +
    theme(legend.text   = element_text(face = "italic"),
          legend.title  = element_text(face = "bold"),
          panel.grid    = element_blank(),
          plot.title    = element_text(face = "bold", size = 12)) +
    labs(title = dist_name,
         x     = paste0("PC1 (", pct[1], "%)"),
         y     = paste0("PC2 (", pct[2], "%)"),
         color = "Species")
}

pcoa_plots <- mapply(
  make_pcoa_plot,
  dist_name = dist_names,
  dist_obj  = dist_list,
  res       = perm_free,
  SIMPLIFY  = FALSE
)

p_combined <- wrap_plots(pcoa_plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("results/D3_upper/plots/24_beta_diversity/PCoA_combined.png",
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
  terms  <- c("Species", "Residual", "Total")

  safe_val <- function(term, col) {
    if (term %in% rownames(res) && col %in% colnames(res)) res[term, col] else NA
  }

  df <- data.frame(
    Term = terms,
    df   = sapply(terms, function(t) safe_val(t, "Df")),
    F    = sapply(terms, function(t) { v <- safe_val(t, "F");  if (is.na(v)) "\u2014" else round(v, 3) }),
    R2   = sapply(terms, function(t) { v <- safe_val(t, "R2"); if (is.na(v)) "\u2014" else round(v, 3) }),
    p    = sapply(terms, function(t) fmt_p_tbl(safe_val(t, "Pr(>F)")))
  )
  colnames(df) <- c("Term", "df", "F", "R\u00b2", "p")

  flextable(df) %>%
    bold(part = "header") %>%
    set_caption(paste0("PERMANOVA (free permutation, sequential SS) \u2014 ", dist_name, " (D3 upper)")) %>%
    footnote(i = 1, j = 5, part = "header",
             value = as_paragraph("Free permutation (999 permutations). No blocking — each individual sampled once."),
             ref_symbols = "a") %>%
    autofit()
}

make_disp_table_sp <- function(perm_sp, dist_name) {
  f <- round(perm_sp$perm$tab["Groups", "F"], 3)
  p <- perm_sp$perm$tab["Groups", "Pr(>F)"]
  df <- data.frame(Group = "Species", F = f, p = fmt_p_tbl(p), stringsAsFactors = FALSE)
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
    body_add_flextable(make_perm_table(perm_free[[d]], d)) %>%
    body_add_par("", style = "Normal") %>%
    body_add_par("PERMDISP (homogeneity of dispersions)", style = "heading 2") %>%
    body_add_flextable(make_disp_table_sp(permdisp_sp[[d]], d)) %>%
    body_add_par("", style = "Normal")
}
print(doc, target = "results/D3_upper/tables/24_beta_diversity/beta_diversity_stats.docx")

cat("\nAll beta diversity outputs saved to:\n")
cat("  results/D3_upper/plots/24_beta_diversity/\n")
cat("  results/D3_upper/tables/24_beta_diversity/\n")
