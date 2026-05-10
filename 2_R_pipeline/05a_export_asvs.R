rm(list = ls())
setwd("D:/Rwd_claude/darevskia_16S_first")
library(dada2)

load("phyloseqs/seqtab_sanger.R")

asv_seqs    <- colnames(seqtab)
asv_ids     <- paste0("ASV", seq_along(asv_seqs))
fasta_lines <- c(rbind(paste0(">", asv_ids), asv_seqs))
writeLines(fasta_lines, "phyloseqs/asvs_sanger.fasta")

cat("Exported", length(asv_seqs), "ASV sequences to phyloseqs/asvs_sanger.fasta\n")
cat("Next: run scripts/1_bioinformatics/05b_vsearch_chimera.sh via WSL\n")
