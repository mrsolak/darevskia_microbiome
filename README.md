# Supplementary Code — *Darevskia* 16S Gut Microbiome Study

This repository contains all bioinformatics and statistical analysis scripts used in the manuscript. Scripts are organised into three folders reflecting the sequential stages of the analysis pipeline.

---

## Study overview

Gut microbiome of two *Darevskia* lizard species (*D. dahli* and *D. portschinskii*) characterised by 16S rRNA amplicon sequencing. Samples were collected from upper and lower intestine compartments of wild-caught individuals. In *D. portschinskii*, both juvenile and adult animals were sampled.

**Sequencing platform:** Element Biosciences AVITI  
**Target region:** 16S rRNA V3–V4  
**Final dataset:** 49 samples, 1,042 ASVs across five analytical subsets (D1–D5)

---

## Software and versions

| Software | Version | Purpose |
|---|---|---|
| Skewer | 0.2.2 | Demultiplexing and adapter trimming |
| DADA2 | 1.x (R package) | ASV inference |
| vsearch | 2.x | Reference-based chimera removal |
| SILVA | release 138.1 | Taxonomic reference database |
| R | 4.5.1 | All statistical analyses |

**Key R packages:** phyloseq, vegan, picante, lme4, emmeans, MaAsLin2, ggplot2, ggVennDiagram, patchwork, flextable, officer

---

## Folder structure

```
1_bioinformatics/    Shell scripts run in Linux (WSL Ubuntu 24.04 or HPC)
2_R_pipeline/        R scripts: ASV inference through phyloseq construction
3_R_analyses/        R scripts: statistical analyses for datasets D1–D5
```

---

## 1. Bioinformatics (`1_bioinformatics/`)

All scripts are bash scripts run under Linux. Input: raw paired-end FASTQ files from the AVITI sequencer.

| Script | Description |
|---|---|
| `02b_skewer_demux_sanger.sh` | Demultiplexes pooled reads by bare primer sequences using Skewer (`-f sanger`). Produces one FASTQ pair per sample. |
| `03b_skewer_trim_sanger.sh` | Trims 16S primer sequences from demultiplexed reads using Skewer (`-f sanger`). |
| `05b_vsearch_chimera.sh` | Removes chimeric ASVs by reference-based detection against the SILVA 138.1 database using vsearch `--uchime_ref`. |

**Note on `-f sanger`:** The quality format flag must be specified explicitly. Auto-detection fails for high-quality AVITI data because uniformly high quality scores produce no characters below ASCII 64, leaving the format ambiguous.

---

## 2. R pipeline (`2_R_pipeline/`)

Run sequentially in R (v4.5.1) on Windows. Each script loads the output of the previous step. All scripts begin with `setwd()` and can be run standalone.

| Script | Input | Output | Description |
|---|---|---|---|
| `00_utils.R` | — | — | Utility functions (`dupl.concensus`, `merge.duplicates`) sourced by downstream scripts |
| `04b_dada2_sanger.R` | Trimmed FASTQs | `seqtab_sanger.R` | DADA2 ASV inference (112 samples → 35,602 ASVs). `multithread = FALSE` required. |
| `05a_export_asvs.R` | `seqtab_sanger.R` | `asvs_sanger.fasta` | Exports ASV sequences to FASTA for chimera checking |
| `06_taxonomy.R` | `non_chimeras_sanger.fasta` | `tax_table_sanger.R` | Assigns taxonomy against SILVA 138.1 using DADA2 `assignTaxonomy` (minBoot = 80) |
| `07_build_phyloseq_sanger.R` | `seqtab_sanger.R`, `tax_table_sanger.R`, `non_chimeras_sanger.fasta` | `PS_raw_sanger.R` | Builds initial phyloseq object; removes ASVs with NA phylum and Chloroplast-order ASVs (112 samples, 14,529 ASVs) |
| `08_merge_metadata_sanger.R` | `PS_raw_sanger.R`, metadata CSV | `PS_meta_sanger.R` | Merges sample metadata; removes sequencing controls (108 samples) |
| `09_merge_duplicates_sanger.R` | `PS_meta_sanger.R` | `PS_dupl_sanger.R`, `PS_merged_sanger.R` | Validates technical duplicates (Procrustes test), merges duplicate pairs by summing reads; applies low-depth filter (≥ 500 reads); saves merged object (49 samples, 1,042 ASVs) |
| `10_metadata_clean_sanger.R` | `PS_merged_sanger.R` | `PS_merged_clean_sanger.R` | Harmonises metadata values (location names, derived columns); final canonical phyloseq |
| `11_first_look_sanger.R` | `PS_merged_clean_sanger.R` | Plots | Exploratory barplots and PCoA |
| `12_build_tree_sanger.R` | `PS_merged_clean_sanger.R` | Phyloseq with tree | Builds a phylogenetic tree from ASV sequences (DECIPHER alignment + phangorn NJ/GTR) and inserts it into the phyloseq object |

---

## 3. Statistical analyses (`3_R_analyses/`)

All scripts load `PS_merged_clean_sanger.R` (the canonical final phyloseq, 49 samples, 1,042 ASVs) and subset it to the relevant dataset at the top of the script. No modified phyloseq objects are saved to disk.

### Analytical datasets

| Dataset | Subset | n samples | n ASVs | Purpose |
|---|---|---|---|---|
| D1 | All samples | 49 | 1,042 | Full dataset — maximum power |
| D2 | Individuals with both upper AND lower intestine | 44 | 1,023 | Strictly paired design |
| D3 | Upper intestine only | 25 | 478 | Species comparison within upper intestine |
| D4 | Lower intestine only | 24 | 815 | Species comparison within lower intestine |
| D5 | *D. portschinskii* only | 28 | 981 | Life stage comparison (Juvenile vs Adult) |

### Script overview

| Script | Dataset | Analysis |
|---|---|---|
| `12b_barplots_D1.R` | D1 | Stacked barplots — phylum-level relative abundance by species and intestine type |
| `13_alpha_diversity_D1.R` | D1 | Alpha diversity (Observed ASVs, Shannon, Faith's PD); LMM ~ Species × sample_type + (1\|individual_ID); emmeans FDR post-hoc |
| `14_beta_diversity_D1.R` | D1 | Beta diversity (Bray-Curtis, Jaccard, Weighted/Unweighted UniFrac); PERMANOVA with strata = individual_ID; PERMDISP |
| `15_venn_diagram_D1.R` | D1 | Venn diagrams of shared/unique ASVs across 4 groups (DD/DP × Upper/Lower) at 3 prevalence thresholds |
| `16_maaslin2_D1.R` | D1 | Differential abundance (MaAsLin2) at genus, family, order, phylum level; ~ Species + sample_type + (1\|individual_ID); reference: *D. dahli* / Lower intestine |
| `17_barplots_D2.R` | D2 | As D1 |
| `18_alpha_diversity_D2.R` | D2 | As D1 (paired design) |
| `19_beta_diversity_D2.R` | D2 | As D1 (strata permutation) |
| `20_venn_diagram_D2.R` | D2 | As D1 |
| `21_maaslin2_D2.R` | D2 | As D1 |
| `22_barplots_D3.R` | D3 | Stacked barplots by species |
| `23_alpha_diversity_D3.R` | D3 | Alpha diversity; LM ~ Species (no random effect); Wilcoxon test |
| `24_beta_diversity_D3.R` | D3 | Beta diversity; PERMANOVA free permutation ~ Species; PERMDISP |
| `25_venn_diagram_D3.R` | D3 | Venn diagram: DD Upper vs DP Upper |
| `26_maaslin2_D3.R` | D3 | Differential abundance; ~ Species only; reference: *D. dahli* |
| `27_barplots_D4.R` | D4 | As D3 |
| `28_alpha_diversity_D4.R` | D4 | As D3 |
| `29_beta_diversity_D4.R` | D4 | As D3 |
| `30_venn_diagram_D4.R` | D4 | Venn diagram: DD Lower vs DP Lower |
| `31_maaslin2_D4.R` | D4 | As D3 |
| `32_barplots_D5.R` | D5 | Stacked barplots by life stage and intestine type |
| `33_alpha_diversity_D5.R` | D5 | Alpha diversity; LMM ~ Life.stage × sample_type + (1\|individual_ID); reference: Adult |
| `34_beta_diversity_D5.R` | D5 | Beta diversity; PERMANOVA strata = individual_ID ~ Life.stage × sample_type; PERMDISP |
| `35_venn_diagram_D5.R` | D5 | Venn diagram: 4 groups (Juvenile/Adult × Upper/Lower) |
| `36_maaslin2_D5.R` | D5 | Differential abundance; ~ Life.stage + sample_type + (1\|individual_ID); reference: Adult / Lower intestine |
| `38_intragroup_variance_D1.R` | D1 | Interindividual variance: betadisper (2-group species centroid), LMM distance ~ Species + sample_type + (1\|individual_ID); pairwise within-species Bray-Curtis (descriptive) |
| `39_intragroup_variance_D2.R` | D2 | As D1 (paired subset) |
| `40_intragroup_variance_D3.R` | D3 | Interindividual variance: betadisper (2-group species centroid), LM distance ~ Species; pairwise within-species Bray-Curtis (descriptive) |
| `41_intragroup_variance_D4.R` | D4 | As D3 (lower intestine) |
| `42_intragroup_variance_D5.R` | D5 | Interindividual variance: betadisper (2-group life-stage centroid), LMM distance ~ Life.stage + sample_type + (1\|individual_ID); pairwise within-life-stage Bray-Curtis (descriptive) |

### Statistical design summary

| Dataset | Alpha model | Beta permutation | Intragroup variance model | MaAsLin2 random effect |
|---|---|---|---|---|
| D1, D2 | LMM, REML=FALSE | Strata = individual_ID | LMM, REML=FALSE | individual_ID |
| D3, D4 | LM | Free (no strata) | LM | None |
| D5 | LMM, REML=FALSE | Strata = individual_ID | LMM, REML=FALSE | individual_ID |

MaAsLin2 settings (all datasets): normalisation = TSS, transform = LOG, min_prevalence = 0.10, min_abundance = 1×10⁻⁴, max_significance = 0.05.

---

## Reproducibility notes

- Random seeds are set with `set.seed(42)` before all permutation tests.
- Colour scheme: *D. dahli* = `#E69F00`, *D. portschinskii* = `#0072B2` (colour-blind friendly).
- D5 life stage colours: Juvenile = `#56B4E9`, Adult = `#D55E00`.
