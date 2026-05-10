rm(list = ls())
setwd("D:/Rwd_claude/darevskia_16S_first")
library(Biostrings)
library(dada2)

fasta <- readDNAStringSet("phyloseqs/non_chimeras_sanger.fasta")

taxa <- assignTaxonomy(as.character(fasta),
                       "raw_files/silva_nr99_v138.1_train_set.fa.gz",
                       multithread = FALSE,
                       minBoot = 80)

save(taxa, file = "phyloseqs/tax_table_sanger.R")
cat("Saved: phyloseqs/tax_table_sanger.R\n")
cat("Dimensions:", dim(taxa), "\n")
