rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(Maaslin2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

sample_data(PS_merged_sanger)$sample_type <- ifelse(
  grepl("upper", sample_data(PS_merged_sanger)$sample_type, ignore.case = TRUE),
  "Upper intestine", "Lower intestine"
)

# D3: upper intestine only — Species only model, no random effects
PS_D3 <- prune_samples(
  sample_data(PS_merged_sanger)$sample_type == "Upper intestine",
  PS_merged_sanger
)
PS_D3 <- prune_taxa(taxa_sums(PS_D3) > 0, PS_D3)
cat(sprintf("D3: %d samples, %d ASVs\n", nsamples(PS_D3), ntaxa(PS_D3)))

MET <- data.frame(as.data.frame(sample_data(PS_D3)), stringsAsFactors = FALSE)

tax_levels <- c("Genus", "Family", "Order", "Phylum")

threshold_configs <- list(
  filtered = list(
    min_prev  = 0.10,
    min_abund = 1e-4,
    label     = "Prevalence \u2265 10%, Abundance \u2265 1\u00d710\u207b\u2074"
  )
)

SIG_Q <- 0.05

dir.create("results/D3_upper/plots/26_maaslin2",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/D3_upper/tables/26_maaslin2", recursive = TRUE, showWarnings = FALSE)

prepare_level <- function(ps, rank) {
  ps_glom   <- tax_glom(ps, taxrank = rank, NArm = TRUE)
  tax_df    <- as.data.frame(tax_table(ps_glom))
  new_names <- as.character(tax_df[[rank]])
  new_names[is.na(new_names)] <- "Unknown"
  taxa_names(ps_glom) <- make.unique(new_names)

  otu_mat <- if (taxa_are_rows(ps_glom)) {
    as.data.frame(t(as.matrix(otu_table(ps_glom))))
  } else {
    as.data.frame(as.matrix(otu_table(ps_glom)))
  }
  otu_mat
}

# Heatmap: single term (Species only)
make_heatmap <- function(res, title, subtitle) {
  sig_feats <- unique(res$feature[res$qval < SIG_Q])

  if (length(sig_feats) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No significant associations (q < ", SIG_Q, ")"),
                 size = 5, hjust = 0.5) +
        theme_void() +
        labs(title = title, subtitle = subtitle) +
        theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
              plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey50"))
    )
  }

  df <- res %>%
    filter(feature %in% sig_feats) %>%
    mutate(
      term  = "Species\n(D. ports. vs D. dahli)",
      stars = case_when(
        qval < 0.001 ~ "***",
        qval < 0.01  ~ "**",
        qval < 0.05  ~ "*",
        TRUE         ~ ""
      )
    ) %>%
    arrange(coef)

  df$feature <- factor(df$feature, levels = unique(df$feature))
  lim <- max(abs(df$coef), na.rm = TRUE)

  ggplot(df, aes(x = term, y = feature, fill = coef)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = stars), size = 3.8, vjust = 0.8, colour = "black") +
    scale_fill_gradient2(
      low = "#0072B2", mid = "white", high = "#D55E00",
      midpoint = 0, limits = c(-lim, lim),
      name = "Coefficient"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw(base_size = 11) +
    theme(
      axis.title        = element_blank(),
      axis.text.x       = element_text(size = 10, face = "bold"),
      axis.text.y       = element_text(size = 9,  face = "italic"),
      panel.grid        = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle     = element_text(size = 9,  hjust = 0.5, colour = "grey50"),
      plot.caption      = element_text(size = 8,  colour = "grey50", hjust = 0)
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      caption  = paste0("Stars: *** q<0.001  ** q<0.01  * q<0.05  . q<", SIG_Q)
    )
}

# Volcano: single term, no facet
make_volcano <- function(res, title, subtitle) {
  df <- res %>%
    mutate(
      direction = case_when(
        qval >= SIG_Q ~ "ns",
        coef  >  0    ~ "positive",
        TRUE          ~ "negative"
      ),
      neg_log_q = -log10(pmax(qval, 1e-10))
    )

  top_labels <- df %>%
    filter(qval < SIG_Q) %>%
    slice_min(qval, n = 15, with_ties = FALSE) %>%
    pull(feature)
  df$label <- ifelse(df$feature %in% top_labels, df$feature, NA_character_)

  q_line <- -log10(SIG_Q)

  ggplot(df, aes(x = coef, y = neg_log_q, colour = direction)) +
    geom_hline(yintercept = q_line, linetype = "dashed", colour = "grey55", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "solid",  colour = "grey80", linewidth = 0.4) +
    geom_point(alpha = 0.75, size = 2) +
    geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                    fontface = "italic", show.legend = FALSE,
                    segment.colour = "grey60", segment.linewidth = 0.3) +
    scale_colour_manual(
      values = c("ns" = "grey75", "positive" = "#D55E00", "negative" = "#0072B2"),
      guide  = "none"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 9,  hjust = 0.5, colour = "grey50"),
      plot.caption     = element_text(size = 8,  colour = "grey50", hjust = 0)
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Coefficient (D. portschinskii vs D. dahli)",
      y        = expression(-log[10](q)),
      caption  = paste0("Dashed line: q = ", SIG_Q,
                        "  |  Reference group: D. dahli",
                        "\n  Orange: enriched in D. portschinskii vs D. dahli (ref)  |  Blue: enriched in D. dahli (ref)")
    )
}

for (rank in tax_levels) {
  cat(sprintf("\n=== Taxonomic level: %s ===\n", rank))

  otu_mat <- prepare_level(PS_D3, rank)
  cat(sprintf("  Features after agglomeration: %d\n", ncol(otu_mat)))

  for (thr_name in names(threshold_configs)) {
    thr <- threshold_configs[[thr_name]]
    cat(sprintf("  Running: %s (%s)\n", thr_name, thr$label))

    tag      <- paste0(tolower(rank), "_", thr_name)
    out_dir  <- file.path("results/D3_upper/tables/26_maaslin2", tag)
    plot_dir <- file.path("results/D3_upper/plots/26_maaslin2",  tag)
    dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    fit <- tryCatch(
      Maaslin2(
        input_data     = otu_mat,
        input_metadata = MET,
        output         = out_dir,
        fixed_effects  = "Species",
        normalization  = "TSS",
        transform      = "LOG",
        min_prevalence = thr$min_prev,
        min_abundance  = thr$min_abund,
        reference        = "Species,Darevskia dahli",
        max_significance = 0.05,
        cores            = 1,
        plot_heatmap     = FALSE,
        plot_scatter     = FALSE
      ),
      error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
    )

    if (is.null(fit)) next

    res_file <- file.path(out_dir, "all_results.tsv")
    if (!file.exists(res_file)) { cat("  No results file found\n"); next }
    res <- read.table(res_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    n_total <- nrow(res)
    n_sig   <- sum(res$qval < SIG_Q, na.rm = TRUE)
    cat(sprintf("  Total associations: %d | Significant (q < %.2f): %d\n",
                n_total, SIG_Q, n_sig))

    fig_title    <- paste0("Differential abundance \u2014 ", rank, " level (D3 upper)")
    fig_subtitle <- thr$label

    p_heat      <- make_heatmap(res, fig_title, fig_subtitle)
    n_sig_feats <- length(unique(res$feature[res$qval < SIG_Q]))
    heat_h      <- min(22, max(4, n_sig_feats * 0.4 + 3))
    ggsave(file.path(plot_dir, "heatmap.png"),
           p_heat, width = 5, height = heat_h, dpi = 300,
           bg = "white", limitsize = FALSE)

    p_volc <- make_volcano(res, fig_title, fig_subtitle)
    ggsave(file.path(plot_dir, "volcano.png"),
           p_volc, width = 8, height = 5.5, dpi = 300, bg = "white")

    cat(sprintf("  Figures saved: %s\n", plot_dir))
  }
}

cat("\n=== All MaAsLin2 D3 analyses complete ===\n")
cat("Results: results/D3_upper/tables/26_maaslin2/\n")
cat("Figures: results/D3_upper/plots/26_maaslin2/\n")
