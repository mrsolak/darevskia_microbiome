rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(picante)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

# D5: D. portschinskii only
PS_D5 <- prune_samples(
  sample_data(PS_merged_sanger)$Species == "Darevskia portschinskii",
  PS_merged_sanger
)
PS_D5 <- prune_taxa(taxa_sums(PS_D5) > 0, PS_D5)
cat(sprintf("D5: %d samples, %d ASVs\n", nsamples(PS_D5), ntaxa(PS_D5)))

dir.create("results/D5_portschinskii/plots/33_alpha_diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D5_portschinskii/tables/33_alpha_diversity", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1: Compute alpha diversity metrics
# ============================================================

alpha_ps <- estimate_richness(PS_D5, measures = c("Shannon", "Observed"))
alpha_ps$sample <- sample_names(PS_D5)

comm_mat <- if (taxa_are_rows(PS_D5)) {
  t(as.matrix(otu_table(PS_D5)))
} else {
  as.matrix(otu_table(PS_D5))
}
pd_res   <- pd(comm_mat, phy_tree(PS_D5), include.root = TRUE)
alpha_ps$FaithPD <- pd_res$PD[match(alpha_ps$sample, rownames(pd_res))]

MET <- data.frame(as.data.frame(sample_data(PS_D5)), stringsAsFactors = FALSE)
MET$sample <- rownames(MET)
MET$Life.stage <- factor(MET$Life.stage, levels = c("Juv", "Adult"))
alpha_df <- left_join(alpha_ps,
                      MET[, c("sample", "Life.stage", "sample_type", "individual_ID")],
                      by = "sample")

write.csv(alpha_df, "results/D5_portschinskii/tables/33_alpha_diversity/alpha_raw.csv",
          row.names = FALSE)

cat("Alpha diversity summary:\n")
print(summary(alpha_df[, c("Shannon", "Observed", "FaithPD")]))

# ============================================================
# Step 2: Linear mixed models
# Life.stage * sample_type fixed; individual_ID random
# ============================================================

metrics <- c("Observed", "Shannon", "FaithPD")

run_lmm <- function(df, metric) {
  fml <- as.formula(paste(metric, "~ Life.stage * sample_type + (1 | individual_ID)"))
  mod <- lmer(fml, data = df, REML = FALSE)
  list(model = mod, anova = anova(mod))
}

lmm_results <- lapply(metrics, function(m) run_lmm(alpha_df, m))
names(lmm_results) <- metrics

cat("\n=== LMM results ===\n")
for (m in metrics) {
  cat("\n---", m, "---\n")
  print(lmm_results[[m]]$anova)
}

# ============================================================
# Step 3: Pairwise emmeans contrasts + combined boxplot
# ============================================================

type_colors   <- c("Upper intestine" = "#009E73", "Lower intestine" = "#CC79A7")
metric_labels <- c("Observed" = "Observed ASVs", "Shannon" = "Shannon index", "FaithPD" = "Faith's PD")

dodge_w <- 0.7
# x = 1: Juv, x = 2: Adult; fill = Upper (left) / Lower (right)
xpos <- c(juv_lower   = 1 - dodge_w/2,
          juv_upper   = 1 + dodge_w/2,
          adult_lower = 2 - dodge_w/2,
          adult_upper = 2 + dodge_w/2)

get_pairwise_p <- function(mod) {
  emm <- emmeans(mod, ~ Life.stage * sample_type)
  pr  <- as.data.frame(pairs(emm, adjust = "fdr"))
  find_p <- function(g1, g2) {
    idx <- which(grepl(g1, pr$contrast, fixed = TRUE) &
                 grepl(g2, pr$contrast, fixed = TRUE))
    if (!length(idx)) return(NA_real_)
    pr$p.value[idx[1]]
  }
  list(
    juv_LU    = find_p("Juv Lower intestine",   "Juv Upper intestine"),
    adult_LU  = find_p("Adult Lower intestine", "Adult Upper intestine"),
    lower_JA  = find_p("Juv Lower intestine",   "Adult Lower intestine"),
    upper_JA  = find_p("Juv Upper intestine",   "Adult Upper intestine")
  )
}

pw_pvals <- lapply(metrics, function(m) get_pairwise_p(lmm_results[[m]]$model))
names(pw_pvals) <- metrics

fmt_stars <- function(p) {
  if (is.na(p)) return("")
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

build_brackets <- function(m, pv) {
  vals  <- alpha_df[[m]]
  y_max <- max(vals, na.rm = TRUE)
  step  <- diff(range(vals, na.rm = TRUE)) * 0.09
  tick  <- step * 0.25
  data.frame(
    metric = metric_labels[m],
    x1    = c(xpos["juv_lower"],   xpos["adult_lower"],  xpos["juv_lower"],   xpos["juv_upper"]),
    x2    = c(xpos["juv_upper"],   xpos["adult_upper"],  xpos["adult_lower"], xpos["adult_upper"]),
    y_br  = c(y_max + step,        y_max + step,         y_max + 2.8*step,    y_max + 4.0*step),
    tick  = tick,
    label = c(fmt_stars(pv$juv_LU), fmt_stars(pv$adult_LU),
              fmt_stars(pv$lower_JA), fmt_stars(pv$upper_JA)),
    stringsAsFactors = FALSE
  )
}

bracket_df        <- do.call(rbind, lapply(metrics, function(m) build_brackets(m, pw_pvals[[m]])))
bracket_df$metric <- factor(bracket_df$metric, levels = metric_labels[metrics])
bracket_df$x_mid  <- (bracket_df$x1 + bracket_df$x2) / 2

cat("\n=== Pairwise emmeans contrasts (FDR-adjusted) ===\n")
for (m in metrics) {
  cat("\n---", m, "---\n")
  pv <- pw_pvals[[m]]
  cat("  Juv Lower vs Upper:    p =", round(pv$juv_LU,   4), "\n")
  cat("  Adult Lower vs Upper:  p =", round(pv$adult_LU, 4), "\n")
  cat("  Lower: Juv vs Adult:   p =", round(pv$lower_JA, 4), "\n")
  cat("  Upper: Juv vs Adult:   p =", round(pv$upper_JA, 4), "\n")
}

fmt_p_compact <- function(p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  p_str <- ifelse(p < 0.001, "p<0.001", paste0("p=", round(p, 3)))
  paste0(p_str, stars)
}

make_strip <- function(metric_label, anova_res) {
  a <- as.data.frame(anova_res)
  paste0(
    metric_label, "\n",
    "Life stage: F=", round(a["Life.stage",       "F value"], 2), " ", fmt_p_compact(a["Life.stage",       "Pr(>F)"]),
    "  |  Sample type: F=", round(a["sample_type",    "F value"], 2), " ", fmt_p_compact(a["sample_type",    "Pr(>F)"]),
    "  |  Interaction: F=", round(a["Life.stage:sample_type", "F value"], 2), " ", fmt_p_compact(a["Life.stage:sample_type", "Pr(>F)"])
  )
}

strip_labels_vec <- setNames(
  sapply(metrics, function(m) make_strip(metric_labels[m], lmm_results[[m]]$anova)),
  metric_labels[metrics]
)

alpha_long <- alpha_df %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(
    metric     = factor(metric, levels = metrics, labels = metric_labels[metrics]),
    Life.stage = factor(Life.stage, levels = c("Juv", "Adult"))
  )

p_combined <- ggplot(alpha_long, aes(x = Life.stage, y = value, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA,
               position = position_dodge(dodge_w), width = 0.6, alpha = 0.85) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge_w),
              size = 1.8, alpha = 0.7, color = "grey30") +
  geom_segment(data = bracket_df,
               aes(x = x1, xend = x2, y = y_br, yend = y_br),
               inherit.aes = FALSE, linewidth = 0.45) +
  geom_segment(data = bracket_df,
               aes(x = x1, xend = x1, y = y_br - tick, yend = y_br),
               inherit.aes = FALSE, linewidth = 0.45) +
  geom_segment(data = bracket_df,
               aes(x = x2, xend = x2, y = y_br - tick, yend = y_br),
               inherit.aes = FALSE, linewidth = 0.45) +
  geom_text(data = bracket_df,
            aes(x = x_mid, y = y_br, label = label),
            inherit.aes = FALSE, vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = type_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.40))) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             labeller = as_labeller(strip_labels_vec)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x        = element_text(size = 11),
        strip.text          = element_text(size = 8.5, lineheight = 1.4,
                                           margin = margin(t = 5, b = 5)),
        panel.grid.major.x  = element_blank(),
        legend.title        = element_text(face = "bold"),
        legend.position     = "right",
        plot.caption        = element_text(size = 7.5, hjust = 0, color = "grey40",
                                           margin = margin(t = 8))) +
  labs(x = "Life stage", y = NULL, fill = "Sample type",
       caption = "Brackets: pairwise post-hoc contrasts from LMM (emmeans, Satterthwaite df, FDR-adjusted within each metric).\n*** p<0.001  ** p<0.01  * p<0.05  ns p\u22650.05")

ggsave("results/D5_portschinskii/plots/33_alpha_diversity/alpha_diversity_combined.png",
       p_combined, width = 8, height = 12, dpi = 300, bg = "white")

# ============================================================
# Step 4: Save LMM results as .docx
# ============================================================

fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
  ifelse(p < 0.01,  "< 0.01",
  ifelse(p < 0.05,  "< 0.05",
         as.character(round(p, 3)))))
}

make_lmm_table <- function(anova_res, metric_name) {
  df_out <- as.data.frame(anova_res)
  df_out$Term <- rownames(df_out)
  df_out <- df_out[, c("Term", "NumDF", "DenDF", "F value", "Pr(>F)")]
  colnames(df_out) <- c("Term", "numDF", "denDF", "F", "p")
  df_out$F <- round(df_out$F, 3)
  df_out$p <- fmt_p(df_out$p)
  df_out$Term <- gsub("Life.stage:sample_type", "Life stage \u00d7 Sample type", df_out$Term)
  df_out$Term <- gsub("Life.stage", "Life stage", df_out$Term)
  df_out$Term <- gsub("sample_type", "Sample type", df_out$Term)
  flextable(df_out) %>%
    bold(part = "header") %>%
    set_caption(paste0("LMM results \u2014 ", metric_name, " (D5 D. portschinskii)")) %>%
    autofit()
}

doc <- read_docx()
for (m in metrics) {
  ft <- make_lmm_table(lmm_results[[m]]$anova, m)
  doc <- doc %>%
    body_add_par(paste("Alpha diversity LMM:", m), style = "heading 2") %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")
}
print(doc, target = "results/D5_portschinskii/tables/33_alpha_diversity/lmm_results.docx")

# ============================================================
# Step 5: Descriptive statistics table
# ============================================================

desc_tbl <- alpha_df %>%
  group_by(Life.stage, sample_type) %>%
  summarise(
    n             = n(),
    Shannon_mean  = round(mean(Shannon), 3),
    Shannon_sd    = round(sd(Shannon), 3),
    Observed_mean = round(mean(Observed), 1),
    Observed_sd   = round(sd(Observed), 1),
    FaithPD_mean  = round(mean(FaithPD), 3),
    FaithPD_sd    = round(sd(FaithPD), 3),
    .groups = "drop"
  )

write.csv(desc_tbl, "results/D5_portschinskii/tables/33_alpha_diversity/descriptive_stats.csv",
          row.names = FALSE)

ft_desc <- flextable(desc_tbl) %>%
  bold(part = "header") %>%
  set_caption("Alpha diversity descriptive statistics by Life stage and Intestine type — D5") %>%
  autofit()

doc2 <- read_docx() %>% body_add_flextable(ft_desc)
print(doc2, target = "results/D5_portschinskii/tables/33_alpha_diversity/descriptive_stats.docx")

cat("\nAll alpha diversity outputs saved to:\n")
cat("  results/D5_portschinskii/plots/33_alpha_diversity/\n")
cat("  results/D5_portschinskii/tables/33_alpha_diversity/\n")
