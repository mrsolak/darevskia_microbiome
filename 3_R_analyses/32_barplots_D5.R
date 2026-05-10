rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

# D5: D. portschinskii only (Life.stage: Juv / Adult)
PS_D5 <- prune_samples(
  sample_data(PS_merged_sanger)$Species == "Darevskia portschinskii",
  PS_merged_sanger
)
PS_D5 <- prune_taxa(taxa_sums(PS_D5) > 0, PS_D5)
cat(sprintf("D5: %d samples, %d ASVs\n", nsamples(PS_D5), ntaxa(PS_D5)))

dir.create("results/D5_portschinskii/plots/32_barplots", recursive = TRUE, showWarnings = FALSE)

make_phylum_df <- function(ps, top_n = 8) {
  ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_phyl <- tax_glom(ps_rel, taxrank = "Phylum")
  mean_abund <- taxa_sums(ps_phyl) / nsamples(ps_phyl)
  n_top <- min(top_n, ntaxa(ps_phyl))
  top_names <- as.vector(tax_table(ps_phyl)[
    names(sort(mean_abund, decreasing = TRUE))[1:n_top], "Phylum"])
  df <- psmelt(ps_phyl)
  df <- df[, c("Sample", "Life.stage", "sample_type", "Phylum", "Abundance")]
  colnames(df)[c(1, 5)] <- c("sample", "rel_abund")
  df$Phylum <- ifelse(df$Phylum %in% top_names, df$Phylum, "Others")
  df %>% group_by(sample, Life.stage, sample_type, Phylum) %>%
    summarise(rel_abund = sum(rel_abund), .groups = "drop")
}

get_phyl_colors <- function(phyla_vec) {
  phyla_no_other <- phyla_vec[phyla_vec != "Others"]
  n <- length(phyla_no_other)
  pal <- if (n <= 12) brewer.pal(max(3, n), "Paired")[1:n]
         else colorRampPalette(brewer.pal(12, "Paired"))(n)
  cols <- setNames(pal, phyla_no_other)
  cols["Others"] <- "#AAAAAA"
  cols
}

df <- make_phylum_df(PS_D5)
phyla_ordered <- df %>%
  group_by(Phylum) %>% summarise(total = sum(rel_abund), .groups = "drop") %>%
  filter(Phylum != "Others") %>% arrange(desc(total)) %>% pull(Phylum)
df$Phylum <- factor(df$Phylum, levels = c(phyla_ordered, "Others"))
phyl_colors <- get_phyl_colors(c(phyla_ordered, "Others"))

# Order samples by Life.stage then sample_type
df$Life.stage <- factor(df$Life.stage, levels = c("Juv", "Adult"))
sample_order <- df %>% select(sample, Life.stage, sample_type) %>% distinct() %>%
  arrange(Life.stage, sample_type, sample) %>% pull(sample)
df$sample <- factor(df$sample, levels = sample_order)

p <- ggplot(df, aes(x = sample, y = rel_abund, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phyl_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  facet_grid(~ Life.stage + sample_type, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "right", legend.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Relative abundance", fill = "Phylum",
       title = sprintf("Phylum-level composition — D5 D. portschinskii (n = %d)", nsamples(PS_D5)))

ggsave("results/D5_portschinskii/plots/32_barplots/barplot_lifestage_sampletype.png",
       p, width = 12, height = 6, dpi = 300, bg = "white")

cat("Barplot saved to results/D5_portschinskii/plots/32_barplots/\n")
cat("Phyla shown:", paste(phyla_ordered, collapse = ", "), "+ Others\n")
