rm(list = ls())

setwd("D:/Rwd_claude/darevskia_16S_first")

library(phyloseq)
library(DECIPHER)
library(Biostrings)
library(ape)
library(phangorn)

load("phyloseqs/PS_merged_clean_sanger.R")  # loads as PS_merged_sanger

# ============================================================
# Step 1: Extract ASV sequences and align via DECIPHER RNA model
# ============================================================
DNA   <- DNAStringSet(refseq(PS_merged_sanger))
RNA   <- RNAStringSet(DNA)

cat("Aligning", length(RNA), "ASV sequences using DECIPHER RNA structural model...\n")
ALIGN <- AlignSeqs(RNA, processors = 1)  # processors=1 for Windows compatibility
ALIGN <- DNAStringSet(ALIGN)

writeXStringSet(ALIGN, "cache/asvs_aligned.fasta")
cat("Aligned FASTA saved to cache/asvs_aligned.fasta\n")

# ============================================================
# Step 2: Build phylogenetic tree with FastTree
# ============================================================
cat("Building tree with FastTree (GTR + gamma model)...\n")

ret <- system2(
  "D:/Rwd_claude/darevskia_16S_first/FastTree.exe",
  args   = c("-nt", "-gtr", "-gamma", "cache/asvs_aligned.fasta"),
  stdout = "phyloseqs/tree_sanger.nwk",
  stderr = "cache/fasttree.log"
)

if (ret != 0) stop("FastTree failed — check cache/fasttree.log")
cat("Tree saved to phyloseqs/tree_sanger.nwk\n")

# ============================================================
# Step 3: Read tree, midpoint-root, add to phyloseq
# ============================================================
tree <- read.tree("phyloseqs/tree_sanger.nwk")
tree <- midpoint(tree)  # midpoint rooting — no outgroup available

cat("Tree summary:\n")
print(tree)

PS_merged_sanger <- merge_phyloseq(PS_merged_sanger, phy_tree(tree))

save(PS_merged_sanger, file = "phyloseqs/PS_merged_clean_sanger.R")
cat("\nPS_merged_clean_sanger.R updated with phylogenetic tree.\n")
cat("Samples:", nsamples(PS_merged_sanger), "\n")
cat("ASVs:   ", ntaxa(PS_merged_sanger), "\n")
cat("Tree tips:", ntaxa(phy_tree(PS_merged_sanger)), "\n")
