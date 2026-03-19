#============================================
# Population Genetics & Phylogenetic Analyses
#                  BY *** SHOROUK ALDEYARBI ***
# Description:
# This script performs:
# - Multiple sequence alignment (COI gene)
# - Genetic distance estimation
# - Distance heatmap visualization
# - Phylogenetic tree reconstruction & annotation
# ============================================

##installaion:------------------------------------------------------
# install.packages(c(
#   "ape",
#   "ggplot2",
#   "ggdendro",
#   "pheatmap",
#   "dplyr"
# ))
# install.packages("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("DECIPHER", "ggtree"))
# 
# install.packages("remotes")     
# remotes::install_github("pievos101/PopGenome")
# Install PopGenome from the GitHub repository
# devtools::install_github("pievos101/PopGenome")
# install.packages("apex")


##**## Loading packages: -----------------------------------------------------
library(ape)
library(DECIPHER)
library(ggplot2)
library(ggdendro)
library(pheatmap)
library(ggtree)
library(dplyr)
library(vegan)
library(pegas)
library(adegenet)
library(hierfstat)
library(seqinr)
library(PopGenome)

# COX DATA loading and alignment   ---------------------------------------------------------------

# Step 1: Load sequences from a FASTA file
input_file <- "DATA/cox.fas"
dna_data <- readDNAStringSet(input_file, format = "fasta")

# Step 2: Convert to DNAStringSet format for compatibility
dna_data <- DNAStringSet(dna_data)

# Step 3: Remove sequences with excessive missing bases (N or -)
valid_sequences <- dna_data[width(dna_data) > 0]  # Ensure sequences have valid length

# Step 4: Perform Multiple Sequence Alignment using ClustalW via DECIPHER
# Ensure sequences are in DNAStringSet format
valid_sequences <- DNAStringSet(valid_sequences)

# Remove gaps ('-') from all sequences
valid_sequences <- gsub("-", "", valid_sequences)

# Convert back to DNAStringSet
valid_sequences <- DNAStringSet(valid_sequences)

# Now perform alignment
aligned_sequences <- AlignSeqs(valid_sequences, processors = 4)

# Save aligned sequences to a FASTA file
writeXStringSet(aligned_sequences, "aligned_sequences.fasta")

# Print alignment summary
print(aligned_sequences)


# Distance heatmap:  ------------------------------------------------------


# Load aligned sequences
aligned <- readDNAStringSet("aligned_sequences.fasta")

# Convert to matrix for distance calculation
aligned_matrix <- as.matrix(aligned)
dist_matrix <- DistanceMatrix(aligned, correction="TN93+F")
# Plot distance heatmap
# Ensure dist_matrix is a valid distance matrix
pheatmap(as.matrix(dist_matrix),
         clustering_distance_rows = "euclidean",  # Explicitly set distance type
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Distance Heatmap (COX1.gene)")

# Adjusted pheatmap code with layout tweaks
pheatmap(as.matrix(dist_matrix),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Distance Heatmap (CO1.gene)",
         legend = TRUE,                # Ensure legend is on
         legend_breaks = pretty(range(dist_matrix)),  # Ensure legend covers all values
         fontsize = 10,
         cellwidth = 20,
         cellheight = 20,
         angle_col = "45",
         border_color = NA,
         filename = "COX_heatmap_output.png",  # Save to file to avoid display cropping
         width = 10,                       # Width in inches
         height = 10)

#  Tree with node label and background colors -----------------------------
  ## get the tree after mega alignmnet and save it as nwk
# Read and root the tree (adjust the outgroup label as needed)
tree <- read.tree("DATA/tree/cox.nwk")
tree_rooted <- root(tree, outgroup = "D.krausi.outgroup", resolve.root = TRUE)

# Extract tip labels with row numbers
p <- ggtree(tree_rooted)
tree_data <- p$data
tip_data <- tree_data %>% filter(isTip)

# Sanitize function to extract clean tips by row range
sanitize_tips <- function(start, end, tree, tip_data) {
  tips <- tip_data$label[start:end]
  tips <- tips[!is.na(tips) & tips %in% tree$tip.label]
  return(tips)
}

# Define clusters by row ranges
node1_tips <- sanitize_tips(1, 6, tree_rooted, tip_data)     # A.amoreuxi (M/A)
node2_tips <- sanitize_tips(7, 9, tree_rooted, tip_data)     # A.amoreuxi (E)
node3_tips <- sanitize_tips(10, 14, tree_rooted, tip_data)   # A.bicolor
node4_tips <- sanitize_tips(15, 16, tree_rooted, tip_data)   # A.bicolor Egy + A.crassicauda Egy
node5_tips <- sanitize_tips(17, 20, tree_rooted, tip_data)   # A.crassicauda India
node6_tips <- sanitize_tips(21, 23, tree_rooted, tip_data)   # Leiurus + Buthus

# Function to get MRCA safely
safe_MRCA <- function(tips, tree) {
  if (length(tips) >= 2) getMRCA(tree, tips) else NA
}
# Get nodes
node1 <- safe_MRCA(node1_tips, tree_rooted)
node2 <- safe_MRCA(node2_tips, tree_rooted)
node3 <- safe_MRCA(node3_tips, tree_rooted)
node4 <- safe_MRCA(node4_tips, tree_rooted)
node5 <- safe_MRCA(node5_tips, tree_rooted)
node6 <- safe_MRCA(node6_tips, tree_rooted)

# Plot with cluster highlights
p <- ggtree(tree_rooted) + 
  geom_tiplab(size = 3)

# Add background colors to clusters
if (!is.na(node1)) p <- p + geom_hilight(node = node1, fill = "#A6CEE3", alpha = 0.3, extend = 0.1)
if (!is.na(node2)) p <- p + geom_hilight(node = node2, fill = "#B2DF8A", alpha = 0.3, extend = 0.1)
if (!is.na(node3)) p <- p + geom_hilight(node = node3, fill = "#FB9A99", alpha = 0.3, extend = 0.1)
if (!is.na(node4)) p <- p + geom_hilight(node = node4, fill = "#FDBF6F", alpha = 0.3, extend = 0.1)
if (!is.na(node5)) p <- p + geom_hilight(node = node5, fill = "#CAB2D6", alpha = 0.3, extend = 0.1)
if (!is.na(node6)) p <- p + geom_hilight(node = node6, fill = "#FFFF99", alpha = 0.3, extend = 0.1)

# Show the final tree
print(p)

ggsave("COXtree_plot.png", plot = p, width = 10, height = 8)
















