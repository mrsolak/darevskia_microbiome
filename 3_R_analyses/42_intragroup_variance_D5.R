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

# Subset: D. portschinskii only
PS_D5 <- subset_samples(PS_merged_sanger, Species == "Darevskia portschinskii")
PS_D5 <- prune_taxa(taxa_sums(PS_D5) > 0, PS_D5)

dir.create("results/D5_portschinskii/plots/42_intragroup_variance", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D5_portschinskii/tables/42_intragroup_variance", recursive = TRUE, showWarnings = FALSE)

MET <- data.frame(as.data.frame(sample_data(PS_D5)), stringsAsFactors = FALSE)
MET$sample_type <- factor(MET$sample_type, levels = c("Upper intestine", "Lower intestine"))
MET$Life.stage  <- factor(MET$Life.stage,  levels = c("Juv", "Adult"))

stage_colors <- c("Juv" = "#56B4E9", "Adult" = "#D55E00")
stage_labels <- c("Juv" = "Juvenile", "Adult" = "Adult")

bc_dist <- phyloseq::distance(PS_D5, method = "bray")

# ===========================================================
# betadisper: 2-group life stage centroid
# ===========================================================

bd <- betadisper(bc_dist, MET$Life.stage)

bd_df <- data.frame(
  distance      = bd$distances,
  Life.stage    = MET[names(bd$distances), "Life.stage"],
  sample_type   = MET[names(bd$distances), "sample_type"],
  individual_ID = MET[names(bd$distances), "individual_ID"]
)
bd_df$Life.stage  <- factor(bd_df$Life.stage,  levels = levels(MET$Life.stage))
bd_df$sample_type <- factor(bd_df$sample_type, levels = levels(MET$sample_type))

# ===========================================================
# LMM: distance ~ Life.stage + sample_type + (1|individual_ID)
# ===========================================================

mod <- lmer(distance ~ Life.stage + sample_type + (1 | individual_ID),
            data = bd_df, REML = FALSE)
mod_sum <- summary(mod)

cat("\n=== LMM: distance to centroid ~ Life.stage + sample_type + (1|individual_ID) ===\n")
print(mod_sum)

coef_df <- as.data.frame(coef(mod_sum))
coef_df <- cbind(Term = rownames(coef_df), coef_df)
rownames(coef_df) <- NULL
names(coef_df) <- c("Term", "Estimate", "Std.Error", "df", "t", "p")
coef_df <- coef_df %>% mutate(across(where(is.numeric), ~ round(.x, 4)))

ft <- flextable(coef_df) %>%
  bold(part = "header") %>%
  autofit() %>%
  set_caption("LMM: distance to life-stage centroid (Bray-Curtis) ~ Life.stage + sample_type + (1|individual_ID). D5 — D. portschinskii 28 samples.")

doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = "results/D5_portschinskii/tables/42_intragroup_variance/lmm_distance_to_centroid_D5.docx")
write.csv(coef_df, "results/D5_portschinskii/tables/42_intragroup_variance/lmm_distance_to_centroid_D5.csv",
          row.names = FALSE)

# ===========================================================
# Plot A: distance to centroid by Life.stage
# ===========================================================

p_stage  <- coef_df$p[grepl("Adult", coef_df$Term)]
p_annot  <- ifelse(p_stage < 0.001, "p < 0.001", paste0("p = ", round(p_stage, 3)))

p_bd <- ggplot(bd_df, aes(x = Life.stage, y = distance, fill = Life.stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
  geom_jitter(aes(color = Life.stage, shape = sample_type), width = 0.12, size = 2.5, alpha = 0.85) +
  scale_fill_manual(values  = stage_colors, labels = stage_labels) +
  scale_color_manual(values = stage_colors, labels = stage_labels) +
  scale_shape_manual(values = c("Upper intestine" = 16, "Lower intestine" = 17)) +
  scale_x_discrete(labels = stage_labels) +
  labs(
    x        = "Life stage",
    y        = "Distance to life-stage centroid (Bray-Curtis)",
    shape    = "Intestine type",
    title    = expression("Interindividual variance — D5 ("*italic("D. portschinskii")*" only)"),
    subtitle = paste0("LMM life-stage effect: ", p_annot,
                      "  |  circles = upper, triangles = lower intestine")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  ) +
  guides(fill = "none", color = "none")

ggsave("results/D5_portschinskii/plots/42_intragroup_variance/distance_to_centroid_D5.png",
       p_bd, width = 6, height = 5.5, dpi = 300, bg = "white")

# ===========================================================
# Plot B: pairwise within-group Bray-Curtis (descriptive)
# ===========================================================

bc_mat <- as.matrix(bc_dist)
snames <- rownames(bc_mat)
group_vec <- MET$Life.stage
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
within_df$Life.stage <- factor(MET[within_df$s1, "Life.stage"], levels = levels(MET$Life.stage))

summary_pw <- within_df %>%
  group_by(Life.stage) %>%
  summarise(n_pairs = n(), mean = mean(distance), sd = sd(distance), median = median(distance),
            .groups = "drop")
write.csv(summary_pw,
          "results/D5_portschinskii/tables/42_intragroup_variance/pairwise_within_summary_D5.csv",
          row.names = FALSE)

p_pw <- ggplot(within_df, aes(x = Life.stage, y = distance, fill = Life.stage)) +
  geom_violin(width = 0.65, alpha = 0.45, trim = TRUE) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.85) +
  scale_fill_manual(values = stage_colors, labels = stage_labels) +
  scale_x_discrete(labels = stage_labels) +
  labs(
    x        = "Life stage",
    y        = "Pairwise Bray-Curtis dissimilarity (within life stage)",
    title    = expression("Interindividual variance — D5 ("*italic("D. portschinskii")*" only)"),
    subtitle = "Descriptive only — formal inference via LMM on distance-to-centroid",
    caption  = "Each value = one pairwise distance between two individuals of the same life stage"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 9, hjust = 0)
  )

ggsave("results/D5_portschinskii/plots/42_intragroup_variance/pairwise_within_D5.png",
       p_pw, width = 6, height = 5.5, dpi = 300, bg = "white")

cat("\nDone — D5 intragroup variance (", nsamples(PS_D5), "samples)\n")
cat("LMM life-stage p =", round(p_stage, 4), "\n")
