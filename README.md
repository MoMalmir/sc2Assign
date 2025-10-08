# sc2Assign v0.1.0
**sc2Assign** is a novel algorithm designed for robust and accurate cell type assignment in single-cell RNA sequencing (scRNA-seq) data. By leveraging weighted affinity scores and cell-cell distance metrics, sc2Assign enhances the identification of known and rare cell types in heterogeneous scRNA-seq datasets. This approach integrates non-marker gene expression data to refine cell type assignment, overcoming the limitations of traditional clustering and marker-gene-based methods.

![](Images/sc2Assign_blockdiagram.png)

This is the official Github repository of the sc2Assign package presented in the Article 

# Installation
For installing our package you first need to install devtoools following the steps below: 
```
install.packages("devtools")
library(devtools)
```
Next you can use the command below to install sc2Assign: 
```
devtools::install_github("momalmir/sc2Assign")
library(sc2Assign)
```

# Data input
You need to provide the algorithm with sparse count matrix, sparse nomalized matrix and the gene marker set

# Output files
The ouput files are sc2Assign and WAffinity prediicted cell type for each single cell

# Usage
```
pred_results <- sc2Assign(count.matrix, norm.matrix, marker.genes, scData, percentile = 0.25, reduction = 'tsne')
```

*in which*

- **count.matrix:** A sparse matrix containing raw gene expression counts. The rows represent genes, and the columns represent individual cells.

- **norm.matrix:**  A sparse matrix of normalized gene expression values, with the same dimensions as the count.matrix.

- **marker.genes:** A list where each element is a vector of marker genes associated with specific cell types. Each entry in the list corresponds to a particular cell type, and the marker genes are used to identify and assign cells to their respective types based on their expression profiles.

- **scData:** A Seurat object that contains cell metadata, including expression data and dimensionality reduction results (e.g., tSNE, UMAP or PCA). This object is used for storing additional information such as cell-type predictions and performing tasks like cell reassignment based on WAffinity and sc2Assign scores.

- **percentile (default = 0.25):** A numeric value between 0 and 1 that specifies the threshold for identifying uncertain cell assignments. Cells with WAffinity scores that fall below this percentile are considered less confident and are subject to reassignment using the sc2Assign method. The default value of 0.25 means that the lowest 25% of cells (in terms of confidence) will be considered for reassignment.

- **reduction (default = 'tsne'):**  A string specifying the dimensionality reduction method to use when reassigning uncertain cells. Options can include 'tsne', 'umap', 'pca', or other methods available in the Seurat object.

# Tutorial
A comprehensive tutorial explaining how to use the **sc2Assign** function for cell type assignment is provided in the **Tutorial** directory of this package. The tutorial includes detailed instructions, example code, and step-by-step guidance to help you understand and apply the **sc2Assign** method to your own datasets.

To access the tutorial, navigate to the **Tutorial** folder in the package directory, or click [here](https://figshare.com/account/articles/27208077?file=58550638) to view the tutorial files.

# Authors
- **Mostafa Malmir** - malmir.edumail@gmail.com
- **Yufang Jin** - yufang.jin@utsa.edu
- **Yidong Chen** - ChenY8@uthscsa.edu

# Cite sc2Assign
