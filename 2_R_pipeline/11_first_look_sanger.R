rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)
library(vegan)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

dir.create("results/D1_all/plots/11_first_look_sanger",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/D1_all/tables/11_first_look_sanger", recursive = TRUE, showWarnings = FALSE)

# --- Simplify sample_type labels ---
sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type), "Upper intestine", "Lower intestine"
)

# --- Sample counts per main variable ---
MET <- as.data.frame(sample_data(PS_merged_sanger))

write.csv(as.data.frame(table(Species = MET$Species, sample_type = MET$sample_type)),
          "results/D1_all/tables/11_first_look_sanger/species_sampletype_counts.csv",
          row.names = FALSE)
write.csv(as.data.frame(table(Species = MET$Species, Life.stage = MET$Life.stage)),
          "results/D1_all/tables/11_first_look_sanger/species_lifestage_counts.csv",
          row.names = FALSE)

cat("=== Sample counts ===\n")
cat("\nSpecies × sample_type:\n"); print(table(MET$Species, MET$sample_type))
cat("\nSpecies × Life.stage:\n");  print(table(MET$Species, MET$Life.stage))

# --- Italic species labels ---
species_labels <- c(
  "Darevskia dahli"         = expression(italic("Darevskia dahli")),
  "Darevskia portschinskii" = expression(italic("Darevskia portschinskii"))
)
species_colors <- c("Darevskia dahli" = "#E69F00", "Darevskia portschinskii" = "#0072B2")

italic_labeller <- labeller(
  Species = c("Darevskia dahli"         = "D. dahli",
              "Darevskia portschinskii" = "D. portschinskii")
)

phyl_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#999999")

# ============================================================
# HELPER: top-5 phyla + Others data frame
# ============================================================
make_phylum_df <- function(ps) {
  ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_phyl <- tax_glom(ps_rel, taxrank = "Phylum")
  mean_abund <- taxa_sums(ps_phyl) / nsamples(ps_phyl)
  top5_names <- as.vector(tax_table(ps_phyl)[
    names(sort(mean_abund, decreasing = TRUE))[1:5], "Phylum"])
  df <- psmelt(ps_phyl)
  df <- df[, c("Sample", "Species", "sample_type", "Life.stage", "Phylum", "Abundance")]
  colnames(df)[c(1, 6)] <- c("sample", "rel_abund")
  df$Phylum <- ifelse(df$Phylum %in% top5_names, df$Phylum, "Others")
  df %>% group_by(sample, Species, sample_type, Life.stage, Phylum) %>%
    summarise(rel_abund = sum(rel_abund), .groups = "drop")
}

# ============================================================
# BARPLOT 1 — Species × sample_type
# ============================================================
df1 <- make_phylum_df(PS_merged_sanger)
sample_order1 <- df1 %>% select(sample, Species, sample_type) %>% distinct() %>%
  arrange(Species, sample_type) %>% pull(sample)
df1$sample <- factor(df1$sample, levels = sample_order1)
df1$Phylum <- factor(df1$Phylum, levels = c(unique(df1$Phylum[df1$Phylum != "Others"]), "Others"))

p1 <- ggplot(df1, aes(x = sample, y = rel_abund, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phyl_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~ Species + sample_type, scales = "free_x", space = "free_x",
             labeller = italic_labeller) +
  theme_bw(base_size = 13) +
  theme(axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        strip.text         = element_text(size = 11, face = "bold"),
        legend.position    = "right",
        legend.title       = element_text(face = "bold"),
        panel.grid.major.x = element_blank()) +
  labs(x = "Sample ID", y = "Relative abundance (%)", fill = "Phylum")

ggsave("results/D1_all/plots/11_first_look_sanger/barplot_species_sampletype.png",
       p1, width = 12, height = 6, dpi = 300, bg = "white")

# ============================================================
# BARPLOT 2 — Species × Life.stage
# ============================================================
df2 <- make_phylum_df(PS_merged_sanger)
sample_order2 <- df2 %>% select(sample, Species, Life.stage) %>% distinct() %>%
  arrange(Species, Life.stage) %>% pull(sample)
df2$sample <- factor(df2$sample, levels = sample_order2)
df2$Phylum <- factor(df2$Phylum, levels = c(unique(df2$Phylum[df2$Phylum != "Others"]), "Others"))

p2 <- ggplot(df2, aes(x = sample, y = rel_abund, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phyl_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~ Species + Life.stage, scales = "free_x", space = "free_x",
             labeller = italic_labeller) +
  theme_bw(base_size = 13) +
  theme(axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        strip.text         = element_text(size = 11, face = "bold"),
        legend.position    = "right",
        legend.title       = element_text(face = "bold"),
        panel.grid.major.x = element_blank()) +
  labs(x = "Sample ID", y = "Relative abundance (%)", fill = "Phylum")

ggsave("results/D1_all/plots/11_first_look_sanger/barplot_species_lifestage.png",
       p2, width = 12, height = 6, dpi = 300, bg = "white")

# ============================================================
# PCoA — Bray-Curtis
# ============================================================
pcoa  <- ordinate(PS_merged_sanger, method = "PCoA", distance = "bray")
evals <- pcoa$values$Eigenvalues
pct   <- round(evals / sum(evals[evals > 0]) * 100, 1)
x_lab <- paste0("PC1 (", pct[1], "%)")
y_lab <- paste0("PC2 (", pct[2], "%)")

df_pcoa <- data.frame(pcoa$vectors[, 1:2],
                      Species     = MET$Species,
                      sample_type = MET$sample_type)
colnames(df_pcoa)[1:2] <- c("PC1", "PC2")

bc_dist      <- phyloseq::distance(PS_merged_sanger, method = "bray")
MET_df       <- data.frame(MET)
perm_species <- adonis2(bc_dist ~ Species,     data = MET_df, permutations = 999)
perm_type    <- adonis2(bc_dist ~ sample_type, data = MET_df, permutations = 999)

fmt_perm <- function(res) {
  f <- round(res$F[1], 2); r2 <- round(res$R2[1], 3); p <- res$`Pr(>F)`[1]
  p_str <- ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, 3)))
  paste0("PERMANOVA\nF = ", f, ", R\u00B2 = ", r2, "\n", p_str)
}

# PCoA by Species
p3 <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3.5) +
  stat_ellipse(aes(group = Species), type = "t", level = 0.95, linewidth = 0.9) +
  scale_color_manual(values = species_colors, labels = species_labels) +
  annotate("text", x = Inf, y = -Inf, label = fmt_perm(perm_species),
           hjust = 1.05, vjust = -0.3, size = 3.5, lineheight = 1.3) +
  theme_bw(base_size = 13) +
  theme(legend.text  = element_text(face = "italic"),
        legend.title = element_text(face = "bold")) +
  labs(x = x_lab, y = y_lab, color = "Species")

ggsave("results/D1_all/plots/11_first_look_sanger/PCoA_by_species.png",
       p3, width = 8, height = 6, dpi = 300, bg = "white")

# PCoA by sample_type
p4 <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3.5) +
  stat_ellipse(aes(group = sample_type), type = "t", level = 0.95, linewidth = 0.9) +
  scale_color_manual(values = c("Upper intestine" = "#009E73", "Lower intestine" = "#CC79A7")) +
  annotate("text", x = Inf, y = -Inf, label = fmt_perm(perm_type),
           hjust = 1.05, vjust = -0.3, size = 3.5, lineheight = 1.3) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_text(face = "bold")) +
  labs(x = x_lab, y = y_lab, color = "Intestine type")

ggsave("results/D1_all/plots/11_first_look_sanger/PCoA_by_sampletype.png",
       p4, width = 8, height = 6, dpi = 300, bg = "white")

cat("\nAll plots saved to results/D1_all/plots/11_first_look_sanger/\n")
cat("\n=== PERMANOVA: Species ===\n"); print(perm_species)
cat("\n=== PERMANOVA: sample_type ===\n"); print(perm_type)
