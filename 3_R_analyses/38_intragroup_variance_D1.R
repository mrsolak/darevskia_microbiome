rm(list = ls())
setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(vegan)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggtext)
library(dplyr)
library(flextable)
library(officer)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

dir.create("results/D1_all/plots/38_intragroup_variance", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D1_all/tables/38_intragroup_variance", recursive = TRUE, showWarnings = FALSE)

MET <- data.frame(as.data.frame(sample_data(PS_merged_sanger)), stringsAsFactors = FALSE)
MET$sample_type <- factor(MET$sample_type, levels = c("Upper intestine", "Lower intestine"))
MET$Species     <- factor(MET$Species,
                           levels = c("Darevskia dahli", "Darevskia portschinskii"))

species_colors <- c("Darevskia dahli" = "#E69F00", "Darevskia portschinskii" = "#0072B2")
species_labels <- c("Darevskia dahli"        = "*D. dahli*",
                    "Darevskia portschinskii" = "*D. portschinskii*")

bc_dist <- phyloseq::distance(PS_merged_sanger, method = "bray")

# ===========================================================
# betadisper: 2-group species centroid
# ===========================================================

bd <- betadisper(bc_dist, MET$Species)

bd_df <- data.frame(
  distance    = bd$distances,
  Species     = MET[names(bd$distances), "Species"],
  sample_type = MET[names(bd$distances), "sample_type"],
  individual_ID = MET[names(bd$distances), "individual_ID"]
)
bd_df$Species     <- factor(bd_df$Species,     levels = levels(MET$Species))
bd_df$sample_type <- factor(bd_df$sample_type, levels = levels(MET$sample_type))

# ===========================================================
# LMM: distance ~ Species + sample_type + (1|individual_ID)
# REML=FALSE for fixed-effect inference
# ===========================================================

mod <- lmer(distance ~ Species + sample_type + (1 | individual_ID),
            data = bd_df, REML = FALSE)
mod_sum <- summary(mod)

cat("\n=== LMM: distance to centroid ~ Species + sample_type + (1|individual_ID) ===\n")
print(mod_sum)

# Save coefficient table as docx
coef_df <- as.data.frame(coef(mod_sum))
coef_df <- cbind(Term = rownames(coef_df), coef_df)
rownames(coef_df) <- NULL
names(coef_df) <- c("Term", "Estimate", "Std.Error", "df", "t", "p")
coef_df <- coef_df %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

ft <- flextable(coef_df) %>%
  bold(part = "header") %>%
  autofit() %>%
  set_caption("LMM: distance to species centroid (Bray-Curtis) ~ Species + sample_type + (1|individual_ID). D1 — all 49 samples.")

doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "results/D1_all/tables/38_intragroup_variance/lmm_distance_to_centroid_D1.docx")
write.csv(coef_df, "results/D1_all/tables/38_intragroup_variance/lmm_distance_to_centroid_D1.csv",
          row.names = FALSE)

# ===========================================================
# Plot A: distance to centroid by Species
# ===========================================================

# Extract p-value for Species term
p_species <- coef_df$p[grepl("portschinskii", coef_df$Term)]
p_annot   <- ifelse(p_species < 0.001, "p < 0.001", paste0("p = ", round(p_species, 3)))

p_bd <- ggplot(bd_df, aes(x = Species, y = distance, fill = Species)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
  geom_jitter(aes(color = Species, shape = sample_type), width = 0.12, size = 2.5, alpha = 0.85) +
  scale_fill_manual(values  = species_colors, labels = species_labels) +
  scale_color_manual(values = species_colors, labels = species_labels) +
  scale_shape_manual(values = c("Upper intestine" = 16, "Lower intestine" = 17)) +
  scale_x_discrete(labels = c("Darevskia dahli"        = "*D. dahli*",
                               "Darevskia portschinskii" = "*D. portschinskii*")) +
  labs(
    x        = NULL,
    y        = "Distance to species centroid (Bray-Curtis)",
    fill     = "Species", color = "Species",
    shape    = "Intestine type",
    title    = paste0("Interindividual variance — D1 (all samples)"),
    subtitle = paste0("LMM species effect: ", p_annot,
                      "  |  circles = upper, triangles = lower intestine")
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.text      = element_markdown(),
    axis.text.x      = element_markdown(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  ) +
  guides(fill = "none", color = "none")

ggsave("results/D1_all/plots/38_intragroup_variance/distance_to_centroid_D1.png",
       p_bd, width = 6, height = 5.5, dpi = 300, bg = "white")

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
          "results/D1_all/tables/38_intragroup_variance/pairwise_within_summary_D1.csv",
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
    title    = "Interindividual variance — D1 (all samples)",
    subtitle = "Descriptive only — formal inference via LMM on distance-to-centroid",
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

ggsave("results/D1_all/plots/38_intragroup_variance/pairwise_within_D1.png",
       p_pw, width = 6, height = 5.5, dpi = 300, bg = "white")

cat("\nDone — D1 intragroup variance\n")
cat("LMM species p =", round(p_species, 4), "\n")
