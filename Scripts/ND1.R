#============================================
# Population Genetics & Phylogenetic Analyses
#                  BY *** SHOROUK ALDEYARBI ***
# Description:
# This script performs:
# - Multiple sequence alignment (NDI gene)
# - Genetic distance estimation
# - Distance heatmap visualization
# - Phylogenetic tree reconstruction & annotation
# ============================================
##**## Loading packages: -----------------------------------------------------

library(ape)
library(DECIPHER)
library(ggplot2)
library(ggdendro)
library(pheatmap)
library(ggtree)
library(dplyr)

# COX DATA loading and alignment   ---------------------------------------------------------------

# Step 1: Load sequences from a FASTA file
input_file <- "DATA/NADH.fas"
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
writeXStringSet(aligned_sequences, "ND1_aligned_sequences.fasta")

# Print alignment summary
print(aligned_sequences)

# Distance heatmap:  ------------------------------------------------------


# Load aligned sequences
aligned <- readDNAStringSet("ND1_aligned_sequences.fasta")

# Convert to matrix for distance calculation
aligned_matrix <- as.matrix(aligned)
dist_matrix <- DistanceMatrix(aligned, correction="TN93+F")
# Plot distance heatmap
# Ensure dist_matrix is a valid distance matrix ( just to visulize it)
pheatmap(as.matrix(dist_matrix),
         clustering_distance_rows = "euclidean",  # Explicitly set distance type
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Distance Heatmap (ND1.gene)")
# Adjusted pheatmap code with layout tweaks
pheatmap(as.matrix(dist_matrix),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Distance Heatmap (ND1.gene)",
         legend = TRUE,                # Ensure legend is on
         legend_breaks = pretty(range(dist_matrix)),  # Ensure legend covers all values
         fontsize = 10,
         cellwidth = 20,
         cellheight = 20,
         angle_col = "45",
         border_color = NA,
         filename = "ND1_heatmap_output.png",  # Save to file to avoid display cropping
         width = 10,                       # Width in inches
         height = 10)

#  Tree with node label and background colors -----------------------------
## get the tree after mega alignmnet and save it as nwk
# Read and root the tree (adjust the outgroup label as needed)
tree <- read.tree("DATA/tree/nd1.nwk")
tree_rooted <- root(tree, outgroup = "D.krausi.", resolve.root = TRUE)

# Extract tip labels with row numbers
p <- ggtree(tree_rooted)
tree_data <- p$data
tip_data <- tree_data %>% filter(isTip)

# Function to get MRCA safely
safe_MRCA <- function(tips, tree) {
  if (length(tips) >= 2) getMRCA(tree, tips) else NA
}
# Define the sanitize_tips function
sanitize_tips <- function(start, end, tree, tip_data) {
  tips <- tip_data$label[start:end]
  tips <- tips[!is.na(tips) & tips %in% tree$tip.label]
  return(tips)
}

# Define clusters by row ranges
node1_tips <- sanitize_tips(1, 3, tree_rooted, tip_data)     # Cluster 1
node2_tips <- sanitize_tips(4, 5, tree_rooted, tip_data)     # Cluster 2
node3_tips <- sanitize_tips(6, 6, tree_rooted, tip_data)     # Cluster 3 (single row)
node4_tips <- sanitize_tips(7, 7, tree_rooted, tip_data)     # Cluster 4 (single row)
node5_tips <- sanitize_tips(8, 9, tree_rooted, tip_data)     # Cluster 5
node6_tips <- sanitize_tips(10, 11, tree_rooted, tip_data)   # Cluster 6

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

# Add background colors to clusters and extend the coverage to better visualize tips
if (!is.na(node1)) p <- p + geom_hilight(node = node1, fill = "#A6CEE3", alpha = 0.3, extend = 0.5)
if (!is.na(node2)) p <- p + geom_hilight(node = node2, fill = "#B2DF8A", alpha = 0.3, extend = 0.5)
if (!is.na(node3)) p <- p + geom_hilight(node = node3, fill = "#FB9A99", alpha = 0.3, extend = 0.5)
if (!is.na(node4)) p <- p + geom_hilight(node = node4, fill = "#FDBF6F", alpha = 0.3, extend = 0.5)
if (!is.na(node5)) p <- p + geom_hilight(node = node5, fill = "#CAB2D6", alpha = 0.3, extend = 0.6)
if (!is.na(node6)) p <- p + geom_hilight(node = node6, fill = "#FFFF99", alpha = 0.3, extend = 0.7)

# Adjust tree width by controlling xlim and adding more space to the margins
p <- p + xlim(c(0, max(tree_rooted$edge.length) * 1.5))  # Extend the tree horizontally
p <- p + theme(legend.position = "none", plot.margin = unit(c(1, 1, 1, 2), "cm"))  # Adjust margins



# Show the final tree
print(p)
ggsave("ND1tree_plot.png", plot = p, width = 10, height = 8)




                                         
                                         