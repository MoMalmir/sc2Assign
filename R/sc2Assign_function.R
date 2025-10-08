#' ---
#' title: sc2Assign package
#' @author M. Malmir, Y. Chen, Y. Jin
#' @seealso See iScience 2024 paper.
#' @name sc2Assign_package
# Load necessary libraries
#' @importFrom stats cor quantile aggregate complete.cases
#' @importFrom Seurat Embeddings
#' @importFrom Matrix Matrix
NULL  # This is necessary to document package-level imports

#' WAffinity_Score Function
#'
#' This function calculates the WAffinity score for each cell based on marker gene expression.
#'
#' @param count.matrix A sparse matrix of gene expression counts (rows = genes, columns = cells).
#' @param count.matrix.markers A sparse matrix of gene expression counts for marker genes.
#' @param weight A numeric vector representing weights for each marker gene.
#'
#' @return A numeric vector of WAffinity scores for each cell.
#' @export
#' @name WAffinity_Score
WAffinity_Score <- function(count.matrix, count.matrix.markers, weight) {
  weight <- unname(weight)
  numGenes <- nrow(count.matrix.markers)
  if (numGenes == 0) {
    stop("WAffinity_Score: No genes found in count.matrix.markers.")
  }
  
  count.matrix.genes.weighted <- count.matrix.markers * weight
  wa_score <- Matrix::colSums(count.matrix.genes.weighted)
  cpmi <- Matrix::colSums(count.matrix)
  wa_score <- log10((wa_score * 1e6) / (cpmi * numGenes) + 1)
  return(wa_score)
}

#' WAffinity_Weight Function
#'
#' This function calculates the WAffinity weights using the Pearson correlation coefficient between marker genes.
#'
#' @param norm.matrix.markers A normalized sparse matrix of marker genes' expression data.
#' @param gene.markers A vector of marker genes.
#'
#' @return A numeric vector of gene scores (weights) for the marker genes.
#' @export
#' @name WAffinity_Weight
WAffinity_Weight <- function(norm.matrix.markers, gene.markers) {
  if (nrow(norm.matrix.markers) == 0) {
    stop("WAffinity_Weight: norm.matrix.markers is empty (no matching genes).")
  }
  if (nrow(norm.matrix.markers) == 1) {
    return(1)
  }
  
  
  # convert sparse to dense for correlation if needed
  if (inherits(norm.matrix.markers, "dgCMatrix")) {
    norm.matrix.markers <- as.matrix(norm.matrix.markers)
  }
  
  correlation_matrix <- cor(t(norm.matrix.markers))
  diag(correlation_matrix) <- 1
  gene_scores <- rowSums(correlation_matrix, na.rm = TRUE) /
    max(rowSums(correlation_matrix, na.rm = TRUE))
  gene_scores[gene_scores < 0] <- 0
  return(gene_scores)
}

#' SC_WAffinity Function
#'
#' This function computes the WAffinity scores for all cell types based on marker gene expression and weighted scores.
#'
#' @param count.matrix A sparse matrix containing gene expression counts (rows = genes, columns = cells).
#' @param norm.matrix A sparse matrix of normalized gene expression values for the same genes and cells.
#' @param marker.genes A list of marker genes for each cell type.
#'
#' @return A list containing WAffinity score table, predictions, and the gene weights used for each cell type.
#' @export
#' @name SC_WAffinity
SC_WAffinity <- function(count.matrix, norm.matrix, marker.genes) {
  wa_score_table <- matrix(nrow = length(marker.genes), ncol = ncol(count.matrix))
  colnames(wa_score_table) <- colnames(count.matrix)
  rownames(wa_score_table) <- names(marker.genes)
  
  weight_list <- vector("list", length(marker.genes))
  names(weight_list) <- names(marker.genes)
  
  for (i in seq_along(marker.genes)) {
    gene.markers <- marker.genes[[i]]
    gene.markers <- gene.markers[gene.markers %in% rownames(count.matrix)]
    
    
    if (length(gene.markers) == 0) {
      next
    }
    
    count.matrix.markers <- count.matrix[gene.markers, , drop = FALSE]
    norm.matrix.markers <- norm.matrix[gene.markers, , drop = FALSE]
    
    if (length(gene.markers) <= 1) {
      weight <- 1
    } else {
      weight <- WAffinity_Weight(norm.matrix.markers, gene.markers)
    }
    
    weight_list[[i]] <- weight
    wa_score_table[i, ] <- WAffinity_Score(count.matrix, count.matrix.markers, weight)
  }
  
  WAffinity_predict <- rownames(wa_score_table)[apply(wa_score_table, 2, which.max)]
  WAffinity_predict <- ifelse(apply(wa_score_table, 2, max) > 0, WAffinity_predict, "unknown")
  
  return(list(score = wa_score_table, predict = WAffinity_predict, weight = weight_list))
}


#' SC_sc2Assign Function
#'
#' This function refines the WAffinity predictions by reassigning cells using the difference between the top two cell type scores.
#'
#' @param WAffinity_tabel A matrix of WAffinity scores for all cells.
#' @param scData A Seurat object containing cell data and dimensionality reduction results.
#' @param percentile The percentile threshold used for identifying uncertain cell assignments.
#' @param reduction A string specifying the dimensionality reduction method (e.g., 'tsne', 'pca').
#'
#' @return A list containing the misclassified cells and their reassigned cell types.
#' @export
#' @name sc2Assign
SC_sc2Assign <- function(WAffinity_tabel, scData, percentile, reduction = 'tsne') {
  get_max_diff_and_rows <- function(column) {
    sorted_indices <- order(column, decreasing = TRUE)
    if (length(sorted_indices) < 2) return(c(NA, NA, NA))
    first_max <- column[sorted_indices[1]]
    second_max <- column[sorted_indices[2]]
    diff <- first_max - second_max
    first_max_row <- rownames(WAffinity_tabel)[sorted_indices[1]]
    second_max_row <- rownames(WAffinity_tabel)[sorted_indices[2]]
    return(c(diff, first_max_row, second_max_row))
  }

  max_diff_and_rows <- apply(WAffinity_tabel, 2, get_max_diff_and_rows)
  max_diff_df <- as.data.frame(t(max_diff_and_rows), stringsAsFactors = FALSE)
  colnames(max_diff_df) <- c("Diff", "First_Max_Row", "Second_Max_Row")

  max_diff_df$Diff <- as.numeric(max_diff_df$Diff)
  max_diff_df <- max_diff_df[complete.cases(max_diff_df), ]

  threshold <- quantile(max_diff_df$Diff, percentile)
  filtered_df <- max_diff_df[max_diff_df$Diff <= threshold, ]
  misclassified_cells <- rownames(filtered_df)

  if (length(misclassified_cells) == 0) {
    stop("No misclassified cells found based on the given percentile threshold.")
  }

  reduc_coords <- Embeddings(scData, reduction = reduction)
  scData@meta.data$reduc_1 <- reduc_coords[,1]
  scData@meta.data$reduc_2 <- reduc_coords[,2]

  confident_cells <- setdiff(rownames(scData@meta.data), misclassified_cells)
  cell_type_means <- aggregate(cbind(reduc_1, reduc_2) ~ WAffinity, data = scData@meta.data[confident_cells, ], FUN = mean)
  existing_types <- cell_type_means$WAffinity

  reassign_cells <- function(cell) {
    cell_coords <- reduc_coords[cell, , drop = FALSE]
    type1 <- filtered_df[cell, "First_Max_Row"]
    type2 <- filtered_df[cell, "Second_Max_Row"]

    if (!(type1 %in% existing_types) || !(type2 %in% existing_types)) {
      return(scData@meta.data[cell, "WAffinity"])
    }

    mean_type1 <- cell_type_means[cell_type_means$WAffinity == type1, c("reduc_1", "reduc_2"), drop = FALSE]
    mean_type2 <- cell_type_means[cell_type_means$WAffinity == type2, c("reduc_1", "reduc_2"), drop = FALSE]

    dist_to_type1 <- sqrt(sum((cell_coords - mean_type1)^2))
    dist_to_type2 <- sqrt(sum((cell_coords - mean_type2)^2))

    if (dist_to_type1 < dist_to_type2) {
      return(type1)
    } else {
      return(type2)
    }
  }

  misclassified <- intersect(misclassified_cells, rownames(scData@meta.data))
  reassigned_types <- sapply(misclassified_cells, reassign_cells)
  names(reassigned_types) <- misclassified_cells

  return(list(misclassified_cells = misclassified_cells, reassigned_types = reassigned_types))
}

#' sc2Assign Function
#'
#' This function performs cell type assignment using WAffinity and sc2Assign methods.
#'
#' @param count.matrix A sparse matrix containing gene expression counts where rows represent genes and columns represent cells.
#' @param norm.matrix A sparse matrix containing normalized gene expression data for the same genes and cells.
#' @param marker.genes A list where each element is a vector of marker genes corresponding to a particular cell type.
#' @param scData A Seurat object containing cell metadata and dimensional reduction results.
#' @param percentile A numeric value specifying the threshold percentile for sc2Assign (default is 0.25).
#' @param reduction A string indicating the dimensionality reduction method to use (default is 'tsne').
#'
#' @return A list containing:
#' \describe{
#'   \item{WAffinity}{A vector of WAffinity-based cell type predictions.}
#'   \item{sc2Assign}{A vector of cell type predictions refined by the sc2Assign method.}
#' }
#' @export
#'
#' @name sc2Assign
sc2Assign <- function(count.matrix, norm.matrix, marker.genes, scData, percentile = 0.25, reduction = 'tsne') {

  # Ensure all marker genes are included in the count matrix
  geneNames <- rownames(count.matrix)
  marker.genes <- lapply(marker.genes, FUN = function(x) {x[x %in% geneNames]})

  # WAffinity computation
  WAffinity_Score_result <- SC_WAffinity(count.matrix, norm.matrix, marker.genes)
  WAffinity_tabel <- WAffinity_Score_result[["score"]]

  # Assign WAffinity prediction to Seurat object
  scData$WAffinity <- WAffinity_Score_result$predict

  # sc2Assign prediction
  sc2Assign_results <- SC_sc2Assign(WAffinity_tabel, scData, percentile, reduction)
  misclassified_cells <- sc2Assign_results$misclassified_cells
  reassigned_types <- sc2Assign_results$reassigned_types

  # Update the reassigned cell types in Seurat object
  scData$sc2Assign <- scData$WAffinity
  for (cell_name in misclassified_cells) {
    scData@meta.data[cell_name, "sc2Assign"] <- reassigned_types[cell_name]
  }

  # Return both WAffinity and sc2Assign predictions
  return(list(WAffinity = scData$WAffinity, sc2Assign = scData$sc2Assign))
}
