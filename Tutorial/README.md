# sc2Assign v0.1.0
**sc2Assign** is a novel algorithm designed for robust and accurate cell type assignment in single-cell RNA sequencing (scRNA-seq) data. By leveraging weighted affinity scores and cell-cell distance metrics, sc2Assign enhances the identification of known and rare cell types in heterogeneous scRNA-seq datasets. This approach integrates non-marker gene expression data to refine cell type assignment, overcoming the limitations of traditional clustering and marker-gene-based methods.

This is the official Github repository of the sc2Assign package presented in the Article 

# Installation
For installing our package you first need to install devtoools following the steps below: 
```
install.packages("devtools")
library(devtools)
```
Next you can use the command below to install sc2Assign: 
```
devtools::install_github("msmalmir/sc2Assign")
library(sc2Assign)
```

# Data input
You need to provide the algorithm with count matrix, nomalized matrix and the gene marker set

# Output files
The ouput files are sc2Assign and WAffinity prediicted cell type for each single cell

# Usage

