rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(picante)
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

# D3: upper intestine only — each individual appears once, no random effect needed
PS_D3 <- prune_samples(
  sample_data(PS_merged_sanger)$sample_type == "Upper intestine",
  PS_merged_sanger
)
PS_D3 <- prune_taxa(taxa_sums(PS_D3) > 0, PS_D3)
cat(sprintf("D3: %d samples, %d ASVs\n", nsamples(PS_D3), ntaxa(PS_D3)))

dir.create("results/D3_upper/plots/23_alpha_diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("results/D3_upper/tables/23_alpha_diversity", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1: Compute alpha diversity metrics
# ============================================================

alpha_ps <- estimate_richness(PS_D3, measures = c("Shannon", "Observed"))
alpha_ps$sample <- sample_names(PS_D3)

comm_mat <- if (taxa_are_rows(PS_D3)) {
  t(as.matrix(otu_table(PS_D3)))
} else {
  as.matrix(otu_table(PS_D3))
}
pd_res   <- pd(comm_mat, phy_tree(PS_D3), include.root = TRUE)
alpha_ps$FaithPD <- pd_res$PD[match(alpha_ps$sample, rownames(pd_res))]

MET <- data.frame(as.data.frame(sample_data(PS_D3)), stringsAsFactors = FALSE)
MET$sample <- rownames(MET)
alpha_df <- left_join(alpha_ps,
                      MET[, c("sample", "Species", "Life.stage", "individual_ID")],
                      by = "sample")

write.csv(alpha_df, "results/D3_upper/tables/23_alpha_diversity/alpha_raw.csv",
          row.names = FALSE)

cat("Alpha diversity summary:\n")
print(summary(alpha_df[, c("Shannon", "Observed", "FaithPD")]))

# ============================================================
# Step 2: Linear model (no random effect — each individual once)
# + Wilcoxon test for pairwise comparison
# ============================================================

metrics <- c("Observed", "Shannon", "FaithPD")

run_lm <- function(df, metric) {
  fml <- as.formula(paste(metric, "~ Species"))
  mod <- lm(fml, data = df)
  anv <- anova(mod)
  wil <- wilcox.test(as.formula(paste(metric, "~ Species")), data = df,
                     exact = FALSE)
  list(model = mod, anova = anv, wilcox_p = wil$p.value)
}

lm_results <- lapply(metrics, function(m) run_lm(alpha_df, m))
names(lm_results) <- metrics

cat("\n=== LM ANOVA results ===\n")
for (m in metrics) {
  cat("\n---", m, "---\n")
  print(lm_results[[m]]$anova)
  cat("Wilcoxon p-value:", round(lm_results[[m]]$wilcox_p, 4), "\n")
}

# ============================================================
# Step 3: Combined boxplot with significance bracket
# ============================================================

species_colors <- c("D. dahli" = "#E69F00", "D. portschinskii" = "#0072B2")
metric_labels  <- c("Observed" = "Observed ASVs", "Shannon" = "Shannon index", "FaithPD" = "Faith's PD")

fmt_stars <- function(p) {
  if (is.na(p)) return("")
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

fmt_p_compact <- function(p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  p_str <- ifelse(p < 0.001, "p<0.001", paste0("p=", round(p, 3)))
  paste0(p_str, stars)
}

# Strip labels: metric name + ANOVA stats
make_strip <- function(metric_label, anova_res, wilcox_p) {
  a <- as.data.frame(anova_res)
  paste0(
    metric_label, "\n",
    "Species: F=", round(a["Species", "F value"], 2),
    " (LM), p=", round(a["Species", "Pr(>F)"], 3),
    "  |  Wilcoxon: ", fmt_p_compact(wilcox_p)
  )
}

strip_labels_vec <- setNames(
  sapply(metrics, function(m) make_strip(metric_labels[m],
                                         lm_results[[m]]$anova,
                                         lm_results[[m]]$wilcox_p)),
  metric_labels[metrics]
)

# Bracket heights — one bracket per metric
build_bracket <- function(m) {
  vals  <- alpha_df[[m]]
  y_max <- max(vals, na.rm = TRUE)
  step  <- diff(range(vals, na.rm = TRUE)) * 0.09
  tick  <- step * 0.25
  p     <- lm_results[[m]]$wilcox_p
  data.frame(
    metric = metric_labels[m],
    x1     = 1,
    x2     = 2,
    y_br   = y_max + step,
    tick   = tick,
    label  = fmt_stars(p),
    stringsAsFactors = FALSE
  )
}

bracket_df        <- do.call(rbind, lapply(metrics, build_bracket))
bracket_df$metric <- factor(bracket_df$metric, levels = metric_labels[metrics])
bracket_df$x_mid  <- (bracket_df$x1 + bracket_df$x2) / 2

alpha_long <- alpha_df %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(
    metric       = factor(metric, levels = metrics, labels = metric_labels[metrics]),
    Species_abbr = factor(sub("Darevskia ", "D. ", Species),
                          levels = c("D. dahli", "D. portschinskii"))
  )

p_combined <- ggplot(alpha_long, aes(x = Species_abbr, y = value, fill = Species_abbr)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.7, color = "grey30") +
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
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             labeller = as_labeller(strip_labels_vec)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                           size = 11, face = "italic"),
        strip.text          = element_text(size = 8.5, lineheight = 1.4,
                                           margin = margin(t = 5, b = 5)),
        panel.grid.major.x  = element_blank(),
        plot.caption        = element_text(size = 7.5, hjust = 0, color = "grey40",
                                           margin = margin(t = 8))) +
  labs(x = NULL, y = NULL,
       caption = "Bracket: Wilcoxon rank-sum test.  Strip: LM ANOVA F-test.\n*** p<0.001  ** p<0.01  * p<0.05  ns p\u22650.05")

ggsave("results/D3_upper/plots/23_alpha_diversity/alpha_diversity_combined.png",
       p_combined, width = 7, height = 12, dpi = 300, bg = "white")

# ============================================================
# Step 4: Save LM results as .docx
# ============================================================

fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
  ifelse(p < 0.01,  "< 0.01",
  ifelse(p < 0.05,  "< 0.05",
         as.character(round(p, 3)))))
}

make_lm_table <- function(anova_res, wilcox_p, metric_name) {
  df_out <- as.data.frame(anova_res)
  df_out$Term <- rownames(df_out)
  df_out <- df_out[, c("Term", "Df", "F value", "Pr(>F)")]
  colnames(df_out) <- c("Term", "df", "F", "p (LM)")
  df_out$F <- round(df_out$F, 3)
  df_out[["p (LM)"]] <- fmt_p(df_out[["p (LM)"]])

  wilcox_row <- data.frame(
    Term = "Wilcoxon rank-sum", df = NA, F = NA_real_,
    `p (Wilcoxon)` = fmt_p(wilcox_p),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(wilcox_row) <- c("Term", "df", "F", "p (LM)")

  df_out <- rbind(df_out, wilcox_row)

  flextable(df_out) %>%
    bold(part = "header") %>%
    set_caption(paste0("LM / Wilcoxon results \u2014 ", metric_name, " (D3 upper intestine)")) %>%
    autofit()
}

doc <- read_docx()
for (m in metrics) {
  ft <- make_lm_table(lm_results[[m]]$anova, lm_results[[m]]$wilcox_p, m)
  doc <- doc %>%
    body_add_par(paste("Alpha diversity:", m), style = "heading 2") %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")
}
print(doc, target = "results/D3_upper/tables/23_alpha_diversity/lm_results.docx")

# ============================================================
# Step 5: Descriptive statistics table
# ============================================================

desc_tbl <- alpha_df %>%
  group_by(Species) %>%
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

write.csv(desc_tbl, "results/D3_upper/tables/23_alpha_diversity/descriptive_stats.csv",
          row.names = FALSE)

ft_desc <- flextable(desc_tbl) %>%
  bold(part = "header") %>%
  set_caption("Alpha diversity descriptive statistics by Species — D3 upper intestine") %>%
  autofit()

doc2 <- read_docx() %>% body_add_flextable(ft_desc)
print(doc2, target = "results/D3_upper/tables/23_alpha_diversity/descriptive_stats.docx")

cat("\nAll alpha diversity outputs saved to:\n")
cat("  results/D3_upper/plots/23_alpha_diversity/\n")
cat("  results/D3_upper/tables/23_alpha_diversity/\n")
