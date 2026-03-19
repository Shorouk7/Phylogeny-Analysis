# ============================================
# Population Genetics & Phylogenetic Analyses
#                  BY *** SHOROUK ALDEYARBI ***
# Description:
# This script performs population genetics,
# distance correlation, and sequence evolution analyses.
# ============================================
## INSTALLATION--------------------------------------------------
#install.packages("phangorn")
# install.packages("pegas")    # For population genetics stats
# install.packages("vegan")     # For Mantel test
# install.packages("ape")       # For genetic distance
# install.packages("seqinr")    # Optional, for sequence handling
#LOADING LIB-----------------------------------------------------
###Loding lib:
library(ape)
library(vegan)
library(pegas)

#*** TEST1: tajima-----------------------------------------------
#Data set: 
# Load your aligned sequences
co1 <- read.dna("C:/Users/sh_de/Desktop/PHYLOGENTIC/PHylogeny paper/Phylogeny/DATA/cox.fas", format = "fasta")
nd1 <- read.dna("C:/Users/sh_de/Desktop/PHYLOGENTIC/PHylogeny paper/Phylogeny/DATA/NADH.fas", format = "fasta")
concat = read.dna("C:/Users/sh_de/Desktop/PHYLOGENTIC/PHylogeny paper/Phylogeny/DATA/CONCAT.fas", format = "fasta")
   #TESTING
tajima.test(aligned_sequences)
tajima.test(nd1)
tajima.test(concat)

#*****************************************************************#--
                          
#*** TEST1:Population Genetic Differentiation (FST Analysis)-----------------------------------------------


 ## the other test:
# install.packages("hierfstat")
# install.packages("adegenet")
# install.packages("hierfstat")

library(hierfstat)
library(adegenet)
library(hierfstat)

# Example: Convert your DNA data into genind object
# Assuming you have a FASTA or STRUCTURE format dataset with population labels

# For example, using existing DNA sequences with population info:
# Let's say "pop_data" is your population grouping vector

# Import FASTA and assign population
library(ape)
dna <- read.dna("aligned_sequences.fasta", format = "fasta")
pop <- factor(c("Pop1", "Pop1", "Pop2", "Pop2"))  # Change based on your dataset

# Convert to genind object
gen <- DNAbin2genind(dna)
gen@pop <- pop

# Convert to hierfstat-compatible format
hf_data <- genind2hierfstat(gen)

# Compute Fst
basic.stats(hf_data)$overall["Fst"]

#*****************************************************************#

#*** TEST2: Mantel Test for Genetic Distance Correlation -----------------------------------------------


# Load the libraries
library(vegan)
library(ape)

# Step 1: Read aligned sequences (replace with your actual file paths)
co1 <- read.dna("DATA/cox.fas", format = "fasta")
nd1 <- read.dna("DATA/NADH.fas", format = "fasta")


rownames(nd1) <- c(
  "L.quinquestriatus.Egy",
  "B.arenicola.Egy",
  "A.bicolor.Egy",
  "A.crassicauda.Egy",
  "A.amoreuxi.Egy",
  "A.crassicauda.IN.1",
  "A.amoreuxi.M.AM1.1",
  "A.amoreuxi.A.AM1",
  "A.amoreuxi.M.AM2.1",
  "A.bicolor.T.B.1",
  "D.krausi.outgroup"
)


# Step 1: Standardize sample names (remove punctuation and make lowercase if needed)
clean_names <- function(names_vec) {
  names_vec <- gsub("\\.+$", "", names_vec)           # Remove trailing dots
  names_vec <- gsub("\\s+", "", names_vec)            # Remove spaces
  names_vec <- gsub("\\.+", ".", names_vec)           # Convert multiple dots to single dot
  return(names_vec)
}


rownames(co1) <- clean_names(rownames(co1))
rownames(nd1) <- clean_names(rownames(nd1))

# Step 2: Find common names between the two datasets
common_names <- intersect(rownames(co1), rownames(nd1))
cat("Number of matching samples:", length(common_names), "\n")
print(common_names)

# Step 3: Subset both alignments to the same samples and order
co1_common <- co1[common_names, ]
nd1_common <- nd1[common_names, ]

co1_common <- co1_common[order(rownames(co1_common)), ]
nd1_common <- nd1_common[order(rownames(nd1_common)), ]

# Step 4: Compute distance matrices
dist_co1 <- dist.dna(co1_common, model = "raw")
dist_nd1 <- dist.dna(nd1_common, model = "raw")

# Step 5: Perform Mantel test
mantel_result <- mantel(dist_co1, dist_nd1, method = "pearson", permutations = 999)
print(mantel_result)


# Convert distance matrices to vectors
co1_vec <- as.vector(dist_co1)
nd1_vec <- as.vector(dist_nd1)

# Create a data frame for plotting
df <- data.frame(CO1 = co1_vec, ND1 = nd1_vec)

# Plot
library(ggplot2)
ggplot(df, aes(x = CO1, y = ND1)) +
  geom_point(color = "darkblue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Correlation Genetic Distances",
    subtitle = "Mantel r = 0.42, p = 0.006",
    x = "CO1 Pairwise Genetic Distance",
    y = "ND1 Pairwise Genetic Distance"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5)
  )

ggsave("CO1_ND1_MantelPlot.png", width = 8, height = 6, dpi = 300)

#*** TEST3: Transition–Transversion Saturation Analysis -----------------------------------------------
library(phangorn)

# Load alignment (e.g., ND1.fasta) as DNAbin
library(ape)
nd1 <- read.dna("DATA/NADH.fas", format = "fasta")

# Compute transition/transversion plot (saturation test)
distances <- dist.dna(nd1)
transitions <- dist.dna(nd1, model = "TS")
transversions <- dist.dna(nd1, model = "TV")

plot(distances, transitions, col="blue", pch=16, main="Transitions vs Distance")
points(distances, transversions, col="red", pch=16)
legend("topleft", legend=c("Transitions", "Transversions"), col=c("blue", "red"), pch=16)



#*** TEST4:Nucleotide Base Composition Analysis -----------------------------------------------


#install.packages("seqinr")
library(seqinr)

nd1_seq <- read.fasta("DATA/NADH.fas")
base_comp <- lapply(nd1_seq, function(x) {
  tab <- table(factor(toupper(x), levels = c("A", "T", "G", "C")))
  tab / sum(tab)
})
base_df <- do.call(rbind, base_comp)
barplot(t(base_df), beside=TRUE, col=c("aquamarine", "deeppink", "royalblue1", "lightyellow1"),
        legend.text=TRUE, main="Base Composition per Taxon",  ylab = "Proportion of Bases")


length(nd1_seq)  # Should return 11

## the Diagram ##

png("base_composition.png", width = 5000, height = 3000, res = 300)

rownames(base_df) <- names(nd1_seq)
par(mar = c(18, 6, 4, 2) + 0.1)

bp <- barplot(t(base_df), beside = TRUE,
              col = c("aquamarine", "deeppink", "royalblue1", "lightyellow4"),
              main = "Base Composition per Taxon",
              cex.main = 2.2,         
              ylab = "Proportion of Bases",
              cex.lab = 1.6,          
              cex.axis = 1.4,         
              names.arg = rownames(base_df),
              las = 2,
              cex.names = 1.6,
              ylim = c(0, 0.9))

# Add custom legend with larger text
legend("topright", legend = c("A", "T", "G", "C"),
       fill = c("aquamarine", "deeppink", "royalblue1", "lightyellow4"),
       cex = 1.8,              
       bty = "n")              

dev.off()
colors()[grep("yellow", colors())]


#*** TEST5: Substitution Saturation Test (Iss Index) -----------------------------------------------


# Load packages
library(seqinr)
library(ggplot2)

# --- Step 1: Read alignment ---
seqs <- read.fasta("DATA/NADH.fas", as.string = FALSE, set.attributes = FALSE)
num_taxa <- length(seqs)
seq_len <- length(seqs[[1]])
cat("Loaded alignment:", num_taxa, "taxa,", seq_len, "sites\n")

# --- Step 2: Shannon entropy function ---
shannon_entropy <- function(column) {
  freqs <- table(column) / length(column)
  -sum(freqs * log(freqs))
}

# --- Step 3: Observed Iss ---
matrix_seqs <- do.call(cbind, seqs)   # convert list -> alignment matrix
entropies <- apply(matrix_seqs, 1, shannon_entropy)
Iss_obs <- mean(entropies) / log(4)   # normalized
cat("Observed Iss:", Iss_obs, "\n")

# --- Step 4: Critical Iss values (approximation) ---
Iss_sym <- 2.0 / num_taxa
Iss_asym <- 4.0 / (num_taxa - 2)
cat("Critical Iss (sym):", Iss_sym, "\n")
cat("Critical Iss (asym):", Iss_asym, "\n")

# --- Step 5: Prepare data for ggplot ---
df <- data.frame(
  Category = c("Observed", "Critical (sym)", "Critical (asym)"),
  Value = c(Iss_obs, Iss_sym, Iss_asym)
)

# --- Step 6: Enhanced plot with ggplot2 ---
ggplot(df, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = round(Value, 3)), vjust = -0.8, size = 4.5) +
  scale_fill_manual(values = c("royalblue", "seagreen", "tomato")) +
  ylim(0, max(df$Value) * 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution Saturation Test",
    y = "Index of Substitution Saturation (Iss)",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# --- Step 6: Enhanced plot with ggplot2 ---
p <- ggplot(df, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = round(Value, 3)), vjust = -0.8, size = 4.5) +
  scale_fill_manual(values = c("royalblue", "seagreen", "tomato")) +
  ylim(0, max(df$Value) * 1.2) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Substitution Saturation Test",
    y = "Index of Substitution Saturation (Iss)",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Show on screen
print(p)

# --- Save in high resolution (example: 600 dpi, 8x6 inches) ---
ggsave("Iss_Test_plot.png", plot = p, width = 8, height = 6, dpi = 600)
ggsave("Iss_Test_plot.pdf", plot = p, width = 8, height = 6)  # vector format (best for journals)






