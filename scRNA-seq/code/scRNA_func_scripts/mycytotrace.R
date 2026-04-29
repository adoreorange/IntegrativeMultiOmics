#cytoTRACE analysis
#using mnn corrected data
library(Seurat)
library(tidyverse)
library(CytoTRACE)
source("/home/adore_org/B_scRNA-seq/analysis/scRNA_scripts/mycytotrace.R")


epi.mnn <- anndata::read_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/epi.mnncorrect.h5ad")
epi.mnn.data <- epi.mnn$X
epi.mnn.data <- t(epi.mnn.data)

epi_cytotracez <- CytoTRACEz(epi.mnn.data)
epi.mnn.data.pheno <- as.character(epi.mnn$obs$Cluster)
names(epi.mnn.data.pheno) <- rownames(epi.mnn$obs)

dir.create("/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/",recursive = T)
plotCytoTRACE(epi_cytotracez, phenotype = epi.mnn.data.pheno,
              outputDir ="/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/")
plotCytoGenes(epi_cytotracez, numOfGenes = 20,
              outputDir = "/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/")
saveRDS(epi_cytotracez,"/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/epi_cytotracez.rds")
#cytoTRACE analysis
#using mnn corrected data
library(Seurat)
library(tidyverse)
library(CytoTRACE2)
CytoTRACEz <- function(mat, 
                       enableFast = TRUE, 
                       ncores = 1, 
                       subsamplesize = 1000){
  cat("input nromalized matrix which has been batch-corrected.\n")  #like mnn output
  cat("this function was formated from CytoTRACE.\n")
  
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  a1 <- mat
  a2 <- NULL #batch=NULL
  batch=NULL
  if (ncol(mat) < 3000) {
    enableFast = FALSE
    message("The number of cells in your dataset is less than 3,000. Fast mode has been disabled.")
  }
  else {
    message("The number of cells in your dataset exceeds 3,000. CytoTRACE will now be run in fast mode (see documentation). You can multi-thread this run using the 'ncores' flag. To disable fast mode, please indicate 'enableFast = FALSE'.")
  }
  
  pqgenes <- is.na(rowSums(mat > 0)) | apply(mat, 1, var) == 0
  num_pqgenes <- length(which(pqgenes == TRUE))
  mat <- mat[!pqgenes, ]
  if (num_pqgenes > 0) {
    warning(paste(num_pqgenes, "genes have zero expression in the matrix and were filtered"))
  }
  
  
  if (enableFast == FALSE) {
    size <- ncol(mat)
  }
  else if (enableFast == TRUE & subsamplesize < ncol(mat)) {
    size <- subsamplesize
  }
  else if (enableFast == TRUE & subsamplesize >= ncol(mat)) {
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }
  
  #sampling
  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat)%%chunk)))
  
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately", 
                round(mean(unlist(lapply(subsamples, length)))), "cells each using", 
                min(chunk, ncores), "/", ncores, "core(s)"))
  message(paste("Pre-processing data and generating similarity matrix..."))
  
  batches <- parallel::mclapply(subsamples,
                                mc.cores = min(chunk,ncores), 
                                function(subsample) {
                                  mat <- mat[, subsample]
                                  batch <- batch[subsample]
                                  
                                  #skpping cell QC: na or nFeature_RNA < 10
                                  #do not normalize data again, because input was normalized data
                                  counts <- apply(mat > 0, 2, sum)
                                  mat2 <- mat
                                  
                                  #select 1000 hvg
                                  mvg <- function(matn) {
                                    A <- matn
                                    n_expr <- rowSums(A > 0) 
                                    
                                    #filter gene express in less than 5% cells
                                    A_filt <- A[n_expr >= 0.05 * ncol(A), ]
                                    
                                    #calculated disp delta/miu, retain only 1000 genes
                                    vars <- apply(A_filt, 1, var)
                                    means <- apply(A_filt, 1, mean)
                                    disp <- vars/means
                                    last_disp <- tail(sort(disp), 1000)[1]
                                    A_filt <- A_filt[disp >= last_disp, ]
                                    
                                    return(A_filt)
                                  }
                                  
                                  mat2.mvg <- mvg(mat2)
                                  rm1 <- colSums(mat2.mvg) == 0
                                  mat2 <- mat2[, !rm1]
                                  counts <- counts[!rm1]
                                  
                                  similarity_matrix_cleaned <- function(similarity_matrix) {
                                    D <- similarity_matrix
                                    cutoff <- mean(as.vector(D))
                                    diag(D) <- 0
                                    D[which(D < 0)] <- 0
                                    D[which(D <= cutoff)] <- 0
                                    Ds <- D
                                    D <- D/rowSums(D)
                                    D[which(rowSums(Ds) == 0), ] <- 0
                                    return(D)
                                  }
                                  
                                  D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))
                                  return(list(mat2 = mat2, counts = counts, D = D))
                                })
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  
  filter <- colnames(a1)[-which(colnames(a1) %in% colnames(mat2))]
  
  if (length(filter) > 0) {
    warning(paste(length(filter), "poor quality cells were filtered based on low or no expression. See 'filteredCells' in returned object for names of filtered cells."))
  }
  
  message("Calculating gene counts signature...")
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x, 
  ], counts))
  names(ds2) <- rownames(mat2)
  
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])), 
  ], 2, mean)
  
  samplesize <- unlist(lapply(lapply(batches, function(x) x$counts), 
                              length))
  
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))
  
  D2 <- lapply(batches, function(x) x$D)
  
  regressed <- function(similarity_matrix_cleaned, score) {
    out <- nnls::nnls(similarity_matrix_cleaned, score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }
  
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9) {
    vals <- score
    v_prev <- rep(vals)
    v_curr <- rep(vals)
    for (i in 1:10000) {
      v_prev <- rep(v_curr)
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% 
                           v_curr) + (1 - ALPHA) * vals
      diff <- mean(abs(v_curr - v_prev))
      if (diff <= 1e-06) {
        break
      }
    }
    return(v_curr)
  }
  
  message("Smoothing values with NNLS regression and diffusion...")
  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, 
                                  function(i) {
                                    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
                                    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
                                    cytotrace <- rank(gcs_diffused)
                                  })
  
  #calculate cytotrace score
  cytotrace <- cytotrace_ranked <- unlist(cytotrace)
  cytotrace <- range01(cytotrace)
  
  cytogenes <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x, ], cytotrace))
  names(cytogenes) <- rownames(mat2)
  
  message("Calculating genes associated with CytoTRACE...")
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2)
  
  cytotrace <- cytotrace[colnames(a1)]
  cytotrace_ranked <- cytotrace_ranked[colnames(a1)]
  gcs <- gcs[colnames(a1)]
  counts <- counts[colnames(a1)]
  mat2 <- t(data.frame(t(mat2))[colnames(a1), ])
  
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2) <- colnames(a1)
  
  message("Done")
  return(list(CytoTRACE = cytotrace, CytoTRACErank = cytotrace_ranked, 
              cytoGenes = sort(cytogenes, decreasing = T), GCS = gcs, 
              gcsGenes = sort(ds2, decreasing = T), Counts = counts, 
              filteredCells = filter, exprMatrix = mat2))
}
