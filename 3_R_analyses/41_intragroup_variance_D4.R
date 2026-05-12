rm(list = ls())
setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggtext)
library(dplyr)
library(flextable)
library(officer)
library(lmerTest)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

# Subset: lower intestine only
PS_D4 <- subset_samples(PS_merged_sanger, sample_type == "Lower intestine")
PS_D4 <- prune_taxa(taxa_sums(PS_D4) > 0, PS_D4)

dir.create("results/D4_lower/plots/41_intragroup_variance", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D4_lower/tables/41_intragroup_variance", recursive = TRUE, showWarnings = FALSE)

MET <- data.frame(as.data.frame(sample_data(PS_D4)), stringsAsFactors = FALSE)
MET$Species <- factor(MET$Species,
                       levels = c("Darevskia dahli", "Darevskia portschinskii"))

species_colors <- c("Darevskia dahli" = "#E69F00", "Darevskia portschinskii" = "#0072B2")
species_labels <- c("Darevskia dahli"        = "*D. dahli*",
                    "Darevskia portschinskii" = "*D. portschinskii*")

bc_dist <- phyloseq::distance(PS_D4, method = "bray")

# ===========================================================
# betadisper: 2-group species centroid
# ===========================================================

bd <- betadisper(bc_dist, MET$Species)

bd_df <- data.frame(
  distance = bd$distances,
  Species  = MET[names(bd$distances), "Species"]
)
bd_df$Species <- factor(bd_df$Species, levels = levels(MET$Species))

# ===========================================================
# LM: distance ~ Species (no random effect — each individual appears once)
# ===========================================================

mod <- lm(distance ~ Species, data = bd_df)
mod_sum <- summary(mod)

cat("\n=== LM: distance to centroid ~ Species (D4 lower intestine) ===\n")
print(mod_sum)

coef_df <- as.data.frame(coef(mod_sum))
coef_df <- cbind(Term = rownames(coef_df), coef_df)
rownames(coef_df) <- NULL
names(coef_df) <- c("Term", "Estimate", "Std.Error", "t", "p")
coef_df <- coef_df %>% mutate(across(where(is.numeric), ~ round(.x, 4)))

p_val   <- coef_df$p[grepl("portschinskii", coef_df$Term)]
p_annot <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))

ft <- flextable(coef_df) %>%
  bold(part = "header") %>%
  autofit() %>%
  set_caption("LM: distance to species centroid (Bray-Curtis) ~ Species. D4 — lower intestine only.")

doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = "results/D4_lower/tables/41_intragroup_variance/lm_distance_to_centroid_D4.docx")
write.csv(coef_df, "results/D4_lower/tables/41_intragroup_variance/lm_distance_to_centroid_D4.csv",
          row.names = FALSE)

summary_bd <- bd_df %>%
  group_by(Species) %>%
  summarise(n = n(), mean = mean(distance), sd = sd(distance), median = median(distance),
            .groups = "drop")
write.csv(summary_bd,
          "results/D4_lower/tables/41_intragroup_variance/betadisper_summary_D4.csv",
          row.names = FALSE)

# ===========================================================
# Plot A: distance to centroid by Species
# ===========================================================

p_bd <- ggplot(bd_df, aes(x = Species, y = distance, fill = Species)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
  geom_jitter(aes(color = Species), width = 0.12, size = 2.5, alpha = 0.85) +
  scale_fill_manual(values  = species_colors, labels = species_labels) +
  scale_color_manual(values = species_colors, labels = species_labels) +
  scale_x_discrete(labels = c("Darevskia dahli"        = "*D. dahli*",
                               "Darevskia portschinskii" = "*D. portschinskii*")) +
  labs(
    x        = NULL,
    y        = "Distance to species centroid (Bray-Curtis)",
    title    = "Interindividual variance — D4 (lower intestine only)",
    subtitle = paste0("LM species effect: ", p_annot)
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_markdown(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  )

ggsave("results/D4_lower/plots/41_intragroup_variance/distance_to_centroid_D4.png",
       p_bd, width = 5.5, height = 5, dpi = 300, bg = "white")

# ===========================================================
# Plot B: pairwise within-group Bray-Curtis (descriptive)
# ===========================================================

bc_mat <- as.matrix(bc_dist)
snames <- rownames(bc_mat)
group_vec <- MET$Species
names(group_vec) <- rownames(MET)

pairs_df <- data.frame(
  s1       = snames[row(bc_mat)[upper.tri(bc_mat)]],
  s2       = snames[col(bc_mat)[upper.tri(bc_mat)]],
  distance = bc_mat[upper.tri(bc_mat)],
  stringsAsFactors = FALSE
)
pairs_df$group1 <- as.character(group_vec[pairs_df$s1])
pairs_df$group2 <- as.character(group_vec[pairs_df$s2])

within_df <- pairs_df[pairs_df$group1 == pairs_df$group2, ]
within_df$Species <- factor(MET[within_df$s1, "Species"], levels = levels(MET$Species))

summary_pw <- within_df %>%
  group_by(Species) %>%
  summarise(n_pairs = n(), mean = mean(distance), sd = sd(distance), median = median(distance),
            .groups = "drop")
write.csv(summary_pw,
          "results/D4_lower/tables/41_intragroup_variance/pairwise_within_summary_D4.csv",
          row.names = FALSE)

p_pw <- ggplot(within_df, aes(x = Species, y = distance, fill = Species)) +
  geom_violin(width = 0.65, alpha = 0.45, trim = TRUE) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.85) +
  scale_fill_manual(values = species_colors, labels = species_labels) +
  scale_x_discrete(labels = c("Darevskia dahli"        = "*D. dahli*",
                               "Darevskia portschinskii" = "*D. portschinskii*")) +
  labs(
    x        = NULL,
    y        = "Pairwise Bray-Curtis dissimilarity (within species)",
    title    = "Interindividual variance — D4 (lower intestine only)",
    subtitle = "Descriptive only — formal inference via Wilcoxon on distance-to-centroid",
    caption  = "Each value = one pairwise distance between two individuals of the same species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_markdown(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 9, hjust = 0)
  )

ggsave("results/D4_lower/plots/41_intragroup_variance/pairwise_within_D4.png",
       p_pw, width = 5.5, height = 5, dpi = 300, bg = "white")

cat("\nDone — D4 intragroup variance (", nsamples(PS_D4), "samples)\n")
cat("LM species p =", round(p_val, 4), "\n")
