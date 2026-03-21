# Scorpion Phylogeny & Population Genetics

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Phylogenetics](https://img.shields.io/badge/Phylogenetics-Analysis-blue?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-Completed-brightgreen?style=for-the-badge)

---

## Author

Assistant Lecturer: Shorouk Aldeyarbi  
Faculty of Science, Port Said University, Egypt  

---

## Related Publication

Evaluating cytochrome C oxidase subunit 1 and NADH dehydrogenase 1 mitochondrial genes for five Buthidae scorpions.  

DOI: 10.7324/JABB.2026.272599  

---

## Overview

This project performs phylogenetic and population genetic analyses of Buthidae scorpions using mitochondrial markers:

- Cytochrome Oxidase I (COI)  
- NADH Dehydrogenase 1 (ND1)  

---

## Methods

- Multiple Sequence Alignment (DECIPHER)  
- Genetic Distance Analysis (TN93 model)  
- Phylogenetic Tree Reconstruction  
- Mantel Test (distance correlation)  
- Substitution Saturation Analysis (Iss)  
- Base Composition Analysis  

---

## Outputs

- Genetic distance heatmaps  
- Phylogenetic trees  
- Mantel correlation plots  
- Saturation test plots  

---

## Requirements

```r
ape
DECIPHER
ggtree
pheatmap
vegan
adegenet
hierfstat
seqinr
```

---

## Project Structure

```
DATA/      → raw sequences and trees  
scripts/   → analysis scripts  
outputs/   → figures and plots  
```

---

## Results

### Phylogenetic Tree (COI)
![COI Tree](Outputs/COXtree_plot.jpg)

---

### Genetic Distance Heatmap
![Heatmap](Outputs/COX_heatmap_output.png)

---

### Mantel Test
![Mantel](Outputs/CO1_ND1_MantelPlot.png)

---

### Substitution Saturation (Iss Test)
![ISS-TEST](Outputs/Iss_Test_plot.png)

---

## Notes

COI demonstrated higher phylogenetic resolution compared to ND1, which showed evidence of substitution saturation.



