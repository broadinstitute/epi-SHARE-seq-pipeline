#DORCs helper function

library(parallel)
library(foreach)
library(chromVAR)
library(Matrix)
library(matrixStats)
library(dplyr)
library(SummarizedExperiment)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(pbmcapply)
library(data.table)
library(tidyft)
library(qlcMatrix)

##custom code from https://github.com/saketkc/blog/blob/main/2022-03-10/SparseSpearmanCorrelation2.ipynb
SparsifiedRanks2 <- function(X) {
  if (class(X)[1] != "dgCMatrix") {
    X <- as(object = X, Class = "dgCMatrix")
  }
  non_zeros_per_col <- diff(x = X@p)
  n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
  offsets <- (n_zeros_per_col - 1) / 2
  x <- X@x
  ## split entries to columns
  col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
  ## calculate sparsified ranks and do shifting
  sparsified_ranks <- unlist(x = lapply(X = seq_along(col_lst), 
                                        FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]))
  ## Create template rank matrix
  X.ranks <- X
  X.ranks@x <- sparsified_ranks
  return(X.ranks)
}

SparseSpearmanCor2 <- function(X, Y = NULL, cov = FALSE) {
  
  # Get sparsified ranks
  rankX <- SparsifiedRanks2(X)
  if (is.null(Y)){
    # Calculate pearson correlation on rank matrices
    return (corSparse(X=rankX, cov=cov))
  }
  rankY <- SparsifiedRanks2(Y)
  return(corSparse( X = rankX, Y = rankY, cov = cov))
}

smoothScoresNN <- function(NNmat,
                           TSSmat,
                           geneList = NULL,
                           barcodesList = NULL,
                           nCores = 1)
{
  if (is.null(rownames(NNmat)))
    stop("NN matrix has to have matching cell IDs as rownames\n")
  if (!all.equal(rownames(NNmat), colnames(TSSmat)))
    stop("Nearest-neighbor matrix and TSS activity score matrix don't have matching cells ..\n")
  cat("Number of cells in supplied TSS matrix: ", ncol(TSSmat),
      "\n")
  cat("Number of genes in supplied TSS matrix: ", nrow(TSSmat),
      "\n")
  cat("Number of nearest neighbors being used per cell for smoothing: ",
      ncol(NNmat),
      "\n")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% rownames(TSSmat)))) {
      cat("One or more of the gene names supplied is not present in the TSS matrix provided: \n")
      cat(geneList[!geneList %in% rownames(TSSmat)], sep = ", ")
      cat("\n")
      stop()
    }
    cat("Running TSS score smoothing for genes:", geneList,
        sep = "\n")
    cat("........\n")
    TSSmat <- TSSmat[rownames(TSSmat) %in% geneList,]
  }
  else {
    if (nrow(TSSmat) > 10000) {
      cat(
        "Running smoothing for all genes in TSS matrix! (n = ",
        nrow(TSSmat),
        ") This is bound to take more time than querying specific markers ..\n",
        sep = ""
      )
    }
  }
  opts <- list()
  pb <- txtProgressBar(min = 0,
                       max = ncol(TSSmat),
                       style = 3)
  progress <- function(n)
    setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  
  
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  
  if (!is.null(barcodesList)) {
    cat("Subsetting to ",
        length(barcodesList),
        " barcodes in dataset..\n")
    NNmat <- NNmat[barcodesList, ]
  }
  cat("Running in parallel using ", nCores, "cores ..\n")
  matL <-
    foreach::foreach(
      x = 1:nrow(NNmat),
      .options.snow = opts,
      .packages = c("Matrix", "data.table", "dplyr")
    ) %dopar% {
      smoothedScore <- data.table(Matrix::rowMeans(TSSmat[, NNmat[x, ]]))
      rownames(smoothedScore) <- rownames(TSSmat)
      colnames(smoothedScore) <- rownames(NNmat)[x]
      smoothedScore
    }
  
  parallel::stopCluster(cl)
  
  close(pb)
  cat("Merging results ..\n")
  smoothedMat <-
    dplyr::bind_cols(matL) %>% data.matrix() %>% Matrix(sparse = TRUE)
  rownames(smoothedMat) <- rownames(TSSmat)
  #stopifnot(all.equal(colnames(smoothedMat), colnames(TSSmat)))
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed),
            "\n"))
  
  return(smoothedMat)
}


extractTFNames <- function(motifIDs) {
  if (all(grepl("_", motifIDs, fixed = TRUE))) {
    sapply(strsplit(sapply(
      strsplit(motifIDs, "_LINE.", fixed = FALSE), "[[", 2
    ), "_", fixed = FALSE), "[[", 2)
  } else {
    message("One or more provided motif IDs do not contain any '_' characters .. returning IDs as is")
    motifIDs
  }
}

all.unique <- function(x) {
  length(x) == length(unique(x))
}


clean_theme <- function() {
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
}

ggcustom <- function(gg,
                     clean = FALSE,
                     splitLeg = TRUE,
                     ...) {
  gg <- gg +  theme(
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.text = element_text(color = "black", size = 5.5),
    axis.title = element_text(size = 7),
    line = element_line(size = 0.235),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    ... # Additional optional parameters passed to theme()
  )
  
  if (clean)
    gg <- gg + clean_theme()
  
  if (splitLeg) {
    leg <- cowplot::get_legend(gg)
    
    gg <- gg + theme(legend.position = "none")
    
    print(gg)
    grid::grid.newpage()
    grid::grid.draw(leg)
  } else {
    print(gg)
  }
  
}


splitAndFetch <- function (vec, delim, part)
{
  if (length(part) == 1) {
    sapply(strsplit(as.character(vec), delim, fixed = TRUE),
           "[[", part)
  }
  else {
    sapply(strsplit(as.character(vec), delim, fixed = TRUE),
           function(x)
             paste(x[part], collapse = delim))
  }
}

centerCounts <- function (obj,
                          doInChunks = TRUE,
                          chunkSize = 1000)
{
  if (!class(obj) %in% c(
    "SummarizedExperiment",
    "RangedSummarizedExperiment",
    "dgCMatrix",
    "dgeMatrix",
    "Matrix"
  ))
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")
  if (ncol(obj) > 10000)
    doInChunks <- TRUE
  if (doInChunks) {
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize,
        " ..\n\n")
    starts <- seq(1, ncol(obj), chunkSize)
  }
  else {
    starts <- 1
  }
  counts.l <- list()
  for (i in 1:length(starts)) {
    beginning <- starts[i]
    if (i == length(starts)) {
      ending <- ncol(obj)
    }
    else {
      ending <- starts[i] + chunkSize - 1
    }
    cat("Computing centered counts for cells: ",
        beginning,
        " to ",
        ending,
        "..\n")
    if (class(obj) == "RangedSummarizedExperiment" | class(obj) ==
        "SummarizedExperiment") {
      m <- SummarizedExperiment::assay(obj[, beginning:ending])
    }
    else {
      m <- obj[, beginning:ending]
    }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    cCounts <- Matrix::t(Matrix::t(m) / cellMeans)
    counts.l[[i]] <- cCounts
    gc()
  }
  cat("Merging results..\n")
  centered.counts <- do.call("cbind", counts.l)
  cat("Done!\n")
  if (class(obj) == "RangedSummarizedExperiment" | class(obj) ==
      "SummarizedExperiment") {
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  }
  else {
    return(centered.counts)
  }
}

chunkCore <- function (A, R, met, pairnames)
{
  cor(x = A, y = R, method = met)[pairnames]
}

fastGenePeakcorr <-
  function (ATAC.se,
            RNAmat,
            genome,
            geneList = NULL,
            windowPadSize = 50000,
            normalizeATACmat = TRUE,
            keepPosCorOnly = TRUE,
            keepMultiMappingPeaks = FALSE,
            nCores = 4,
            n_bg = 100,
            n_BgPairs = 1e+05,
            n_bootstrapPairs = 1e+05,
            p.cut = NULL,
            chunkSize = 50000)
  {
    require(dplyr)
    stopifnot(inherits(ATAC.se, "RangedSummarizedExperiment"))
    stopifnot(inherits(RNAmat, c("Matrix", "matrix")))
    if (!all.equal(ncol(ATAC.se), ncol(RNAmat)))
      stop("Input ATAC and RNA objects must have same number of cells")
    rownames(ATAC.se) <- paste0("Peak", 1:nrow(ATAC.se))
    peakRanges.OG <- granges(ATAC.se)
    if (is.null(rownames(RNAmat)))
      stop("RNA matrix must have gene names as rownames")
    if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
      message("Peaks with 0 accessibility across cells exist ..")
      message("Removing these peaks prior to running correlations ..")
      peaksToKeep <- Matrix::rowSums(assay(ATAC.se)) != 0
      ATAC.se <- ATAC.se[peaksToKeep, ]
      message("Important: peak indices in returned gene-peak maps are relative to original input SE")
    }
    ATACmat <- assay(ATAC.se)
    if (normalizeATACmat)
      ATACmat <- centerCounts(ATACmat, chunkSize = chunkSize)
    peakRanges <- granges(ATAC.se)
    if (any(Matrix::rowSums(RNAmat) == 0)) {
      message("Genes with 0 expression across cells exist ..")
      message("Removing these genes prior to running correlations ..")
      genesToKeep <- Matrix::rowSums(RNAmat) != 0
      RNAmat <- RNAmat[genesToKeep, ]
    }
    cat("Number of peaks in ATAC data:", nrow(ATACmat), "\n")
    cat("Number of genes in RNA data:", nrow(RNAmat), "\n")
    if (!genome %in% c("hg19", "hg38", "mm10"))
      stop(
        "You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n"
      )
    switch(genome,
           hg19 = {
             TSSg <- hg19TSSRanges
           },
           hg38 = {
             TSSg <- hg38TSSRanges
           },
           mm10 = {
             TSSg <- mm10TSSRanges
           })
    names(TSSg) <- as.character(TSSg$gene_name)
    if (!is.null(geneList)) {
      if (length(geneList) == 1)
        stop("Please specify more than 1 valid gene symbol")
      if (any(!geneList %in% names(TSSg))) {
        cat(
          "One or more of the gene names supplied is not present in the TSS annotation specified: \n"
        )
        cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
        cat("\n")
        stop()
      }
      TSSg <- TSSg[geneList]
    }
    genesToKeep <- intersect(names(TSSg), rownames(RNAmat))
    cat(
      "\nNum genes overlapping TSS annotation and RNA matrix being considered: ",
      length(genesToKeep),
      "\n"
    )
    RNAmat <- RNAmat[genesToKeep, ]
    TSSg <- TSSg[genesToKeep]
    TSSflank <- GenomicRanges::flank(TSSg, width = windowPadSize,
                                     both = TRUE)
    cat("\nTaking peak summits from peak windows ..\n")
    peakSummits <- resize(peakRanges, width = 1, fix = "center")
    cat("Finding overlapping peak-gene pairs ..\n")
    genePeakOv <-
      findOverlaps(query = TSSflank, subject = peakSummits)
    numPairs <- length(genePeakOv)
    cat("Found ",
        numPairs,
        "total gene-peak pairs for given TSS window ..\n")
    cat("Number of peak summits that overlap any gene TSS window: ",
        length(unique(subjectHits(genePeakOv))),
        "\n")
    cat("Number of gene TSS windows that overlap any peak summit: ",
        length(unique(queryHits(genePeakOv))),
        "\n\n")
    cat("Determining background gene-peak pairs ..\n")
    if (is.null(rowData(ATAC.se)$bias)) {
      if (genome %in% "hg19")
        myGenome <-
          BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
      if (genome %in% "mm10")
        myGenome <-
          BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
      if (genome %in% "hg38")
        myGenome <-
          BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
      ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
    }
    cat("Using ", n_BgPairs, " background gene-peak pairs ..\n\n")
    cat("Finding ", n_bg, " nearest neighbors ..\n\n")
    
    start.time <- Sys.time()
    set.seed(123)
    bgPairGenes <- sample(1:length(genesToKeep), n_BgPairs,
                          replace = TRUE)
    bgPairPeaks <- sample(1:length(peakSummits), n_BgPairs,
                          replace = TRUE)
    bgPairFeatures <-
      data.table(
        GC = ATAC.se@rowRanges$bias[bgPairPeaks],
        accessibility = Matrix::rowMeans(ATACmat)[bgPairPeaks],
        expression = Matrix::rowMeans(RNAmat)[bgPairGenes]
      )
    obPairGenes <- genePeakOv@from
    obPairPeaks <- genePeakOv@to
    obPairFeatures <-
      data.table(
        GC = ATAC.se@rowRanges$bias[obPairPeaks],
        accessibility = Matrix::rowMeans(ATACmat)[obPairPeaks],
        expression = Matrix::rowMeans(RNAmat)[obPairGenes]
      )
    allPairFeatures <- scale(rbind(as.matrix(bgPairFeatures),as.matrix(obPairFeatures)))
    bgPairFeatures <- allPairFeatures[1:nrow(bgPairFeatures), ]
    obPairFeatures <- allPairFeatures[(nrow(bgPairFeatures) + 1):(nrow(bgPairFeatures) + nrow(obPairFeatures)), ]
    bgPairsInds <-
      FNN::get.knnx(
        data = as.data.frame(bgPairFeatures),
        query = as.data.frame(obPairFeatures),
        k = n_bg
      )$nn.index
    
    uniq_ids = unique(unlist(split(bgPairsInds, seq(nrow(bgPairsInds))),recursive=F))
    
    end.time <- Sys.time()
    
    cat(paste(
      "\nTime Elapsed to find nearest neighbors: ",
      round(end.time - start.time, 2),
      units(end.time - start.time),
      "\n"
    ))
    
    cat(paste(
      "\n# of unique nearest neighbors: ",
      length(uniq_ids),
      "\n"
    ))
    
    library(doParallel)
    
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    ##Calculate corr for bg pairs
    corPairs <- data.table(Gene = bgPairGenes[uniq_ids],
                           Peak = bgPairPeaks[uniq_ids],
                           ID = uniq_ids,
                           stringsAsFactors = FALSE)
    
    #corPairs_grp = corPairs[, .(Peak = .(Peak)), by = .(Gene)]
    #numPairs <- nrow(corPairs_grp)
    #pairsPerChunk = as.integer(numPairs/100)
   
    pairsPerChunk = 1000
    numPairs <- length(uniq_ids)
    
    starts <- seq(1, numPairs, pairsPerChunk)
    ends <- starts + pairsPerChunk - 1
    ends[length(ends)] <- numPairs
    
    chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
    
    start.time <- Sys.time()
    
    corList_bg <-
      foreach(
        x = 1:length(chunkList),
        .combine = c,
        .inorder = T,
        .multicombine = TRUE,
        .maxcombine = 1000,
        .export = c("SparseSpearmanCor2","SparsifiedRanks2"),
        .packages = c("dplyr",
                      "Matrix","qlcMatrix")
      ) %dopar% {
        chunk = chunkList[[x]]
        geneIndices <- corPairs$Gene[chunk[1]:chunk[2]]
        peakIndices <- corPairs$Peak[chunk[1]:chunk[2]]
        
        diag(SparseSpearmanCor2(Matrix::t(ATACmat[peakIndices, , drop = FALSE]), Matrix::t(RNAmat[geneIndices, , drop = FALSE]))) 
      }
    
    end.time <- Sys.time()
    
    cat(paste(
      "\nTime Elapsed to calculate background corr: ",
      round(end.time - start.time, 2),
      units(end.time - start.time),
      "\n"
    ))

    parallel::stopCluster(cl)
    
    corList_bg = data.table(ID = uniq_ids, corr = corList_bg) #maintain background IDs
    setkey(corList_bg,ID) #assign IDs as keys for faster search/subset
    
    ##Calculate corr for ob pairs
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
     corPairs <- data.table(Gene = obPairGenes,
                            Peak = obPairPeaks,
                            stringsAsFactors = FALSE)
     
     #corPairs_grp = corPairs[, .(Peak = .(Peak)), by = .(Gene)]
     #numPairs <- nrow(corPairs_grp)
     #pairsPerChunk = as.integer(numPairs/100)
     
     pairsPerChunk = 1000
     numPairs <- length(obPairGenes)
     
     starts <- seq(1, numPairs, pairsPerChunk)
     ends <- starts + pairsPerChunk - 1
     ends[length(ends)] <- numPairs
     
     chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
     
     start.time <- Sys.time()
     
     # corList_ob <- foreach(
     #   x = 1:nrow(corPairs_grp),
     #   .combine = c,
     #   .inorder = T,
     #   .multicombine = TRUE,
     #   .export = c("SparseSpearmanCor2","SparsifiedRanks2"),
     #   .packages = c("dplyr",
     #                 "Matrix","qlcMatrix")
     # ) %dopar% {
     #   tryCatch({
     #     geneIndices <- corPairs_grp[x,]$Gene
     #     peakIndices <- corPairs_grp[x,]$Peak[[1]]
     #     
     #     sub_mat = cbind(Matrix::t(RNAmat[geneIndices, , drop = FALSE]), Matrix::t(ATACmat[peakIndices, , drop = FALSE]))
     #     
     #     sapply(2:ncol(sub_mat), function(i) SparseSpearmanCor2(sub_mat[,1, drop=F], sub_mat[,i, drop=F]))
     #     
     #   },
     #   error = function(e){ 
     #     print(e)
     #   })
     # }
     
     corList_ob <-
       foreach(
         x = 1:length(chunkList),
         .combine = c,
         .inorder = T,
         .multicombine = TRUE,
         .maxcombine = 1000,
         .export = c("SparseSpearmanCor2","SparsifiedRanks2"),
         .packages = c("dplyr",
                       "Matrix","qlcMatrix")
       ) %dopar% {
         chunk = chunkList[[x]]
         geneIndices <- corPairs$Gene[chunk[1]:chunk[2]]
         peakIndices <- corPairs$Peak[chunk[1]:chunk[2]]
         
         diag(SparseSpearmanCor2(Matrix::t(ATACmat[peakIndices, , drop = FALSE]), Matrix::t(RNAmat[geneIndices, , drop = FALSE]))) 
       }
    
     end.time <- Sys.time()
     
     cat(paste(
       "\nTime Elapsed to calculate observed corr: ",
       round(end.time - start.time, 2),
       units(end.time - start.time),
       "\n"
     ))
     
     parallel::stopCluster(cl)
    
    start.time <- Sys.time()
    
    pvals <- pbmcapply::pbmcmapply(function(pair_ind) {
      obcor <- corList_ob[pair_ind]
      bgCorrs <- corList_bg[ID %in% bgPairsInds[pair_ind, ]]$corr
      pval <- 1 - stats::pnorm(q = obcor,
                               mean = mean(bgCorrs),
                               sd = sd(bgCorrs))
    }, 1:length(obPairGenes), mc.cores = nCores)
    time_elapsed <- Sys.time() - start.time
    cat(paste(
      "\nTime Elapsed to calculate pvals: ",
      round(time_elapsed, 2),
      units(time_elapsed),
      "\n"
    ))
    corrResults <- data.table(
      Peak = obPairPeaks,
      Gene = obPairGenes,
      rObs = corList_ob,
      pvalZ = pvals,
      stringsAsFactors = FALSE
    )
    if (keepPosCorOnly) {
      cat("Only considering positive associations ..\n")
      #corrResults <- corrResults %>% filter(rObs > 0)
      corrResults = corrResults[rObs > 0]
      
    }
    if (!keepMultiMappingPeaks) {
      cat("Keeping max correlation for multi-mapping peaks ..\n")
      #corrResults <-
      # corrResults %>% group_by(Peak) %>% filter(rObs ==
      #                                            max(rObs))
      corrResults = corrResults[, .SD[rObs == max(rObs)], by = .(Peak)]
    }
    corrResults$Gene <-
      as.character(TSSg$gene_name)[corrResults$Gene]
    corrResults$Peak <-
      as.numeric(splitAndFetch(rownames(ATACmat)[corrResults$Peak],
                               "Peak", 2))
    cat("\nFinished!\n")
    if (!is.null(p.cut)) {
      cat("Using significance cut-off of ",
          p.cut,
          " to subset to resulting associations\n")
      corrResults <- corrResults[pvalZ <= p.cut, ]
    }
    corrResults$PeakRanges <-
      paste(as.character(seqnames(peakRanges.OG[corrResults$Peak])),
            paste(start(peakRanges.OG[corrResults$Peak]), end(peakRanges.OG[corrResults$Peak]),
                  sep = "-"), sep = ":")
    #return(corrResults %>% as.data.frame(stringsAsFactors = FALSE) %>%
    #        dplyr::select(c(
    #         "Peak", "PeakRanges", "Gene", "rObs",
    #        "pvalZ"
    #     )))
    return(as.data.frame(corrResults[, .(Peak, PeakRanges, Gene, rObs, pvalZ)]))
  }

# Function to make J plot of significant peak-gene assocoations to call DORCs using
dorcJPlot <-
  function(dorcTab,
           # table returned from runGenePeakcorr function
           cutoff = 7,
           labelTop = 25,
           returnGeneList = FALSE, # Returns genes passing numPeak filter
           cleanLabels = TRUE,
           labelSize = 4,
           ... ) { # Additional params passed to ggrepel
    stopifnot(all(c("Peak", "Gene", "pvalZ") %in% colnames(dorcTab)))
    
    # Count the number of significant peak associations for each gene (without pre-filtering genes)
    numDorcs <-
      dorcTab  %>% dplyr::group_by(Gene) %>% dplyr::tally() %>% dplyr::arrange(desc(n))
    numDorcs$Index <- 1:nrow(numDorcs) # Add order index
    numDorcs %>% as.data.frame(stringsAsFactors = FALSE) -> numDorcs
    rownames(numDorcs) <- numDorcs$Gene
    
    dorcGenes <- numDorcs$Gene[numDorcs$n >= cutoff]
    
    numDorcs <- numDorcs %>%
      dplyr::mutate(isDORC = ifelse(Gene %in% dorcGenes, "Yes", "No")) %>%
      dplyr::mutate(Label = ifelse(Gene %in% dorcGenes[1:labelTop], Gene, ""))
    
    # Plot
    dorcG <-
      ggplot(numDorcs, aes(
        x = Index,
        y = n,
        color = isDORC,
        label = Label
      )) +
      geom_hline(linetype = "dotted", yintercept = cutoff) +
      geom_vline(linetype = "dotted", xintercept = max(numDorcs[numDorcs$Gene %in% dorcGenes, "Index"])) +
      geom_point(size = 0.8) +
      geom_line() +
      scale_color_manual(values = c("gray65", "firebrick")) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      theme_classic() +
      labs(
        y = "Number of correlated peaks",
        x = "Ranked genes",
        title = paste0("# DORCs: ( n >= ", cutoff, ") = ", length(dorcGenes))
      ) +
      theme(
        axis.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_x_reverse() # flip so we can add labels later, if needed, with more space
    
    if (cleanLabels) {
      dorcG <-
        dorcG + ggrepel::geom_label_repel(
          size = labelSize,
          max.iter = 100,
          max.overlaps = Inf,
          fontface = "italic",
          ...
        )
    } else {
      dorcG <-
        dorcG + ggplot2::geom_text(size = labelSize, fontface = "italic", ...)
    }
    print(dorcG)
    if (returnGeneList)
      return(numDorcs)
    else
      return(dorcG)
  }


getCountsFromFrags <- function (fragFile,
                                peakFile,
                                barcodeList = NULL,
                                maxFragLength = NULL,
                                addColData = TRUE) {
  
  myFrags = fread(fragFile, sep="\t", header=F, data.table=FALSE)
  myPeaks = fread(peakFile, sep="\t", header=F, data.table=FALSE)
  
  peaks = makeGRangesFromDataFrame(myPeaks,seqnames.field = "V1",start.field = "V2",end.field = "V3",starts.in.df.are.0based = TRUE)
  GA = makeGRangesFromDataFrame(myFrags, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  rm(myFrags)
  rm(myPeaks)
  start_time <- Sys.time()
  
  #GA <- fragRanges
  if (!is.null(barcodeList)) {
    cat("Retaining only select barcodes specified within list ..\\n")
    GA <- GA[mcols(GA)[, 1] %in% barcodeList]
  }
  
  if (ncol(mcols(GA)) > 1) {
    if (inherits(mcols(GA)[2][, 1], "character")) {
      colnames(mcols(GA)) <- c("barcodeID", "readID")
    }
    else if (inherits(mcols(GA)[2][, 1], "integer")) {
      colnames(mcols(GA)) <- c("barcodeID", "pcrDup")
    }
  }
  else {
    colnames(mcols(GA)) <- "barcodeID"
  }
  if (!is.null(maxFragLength)) {
    cat("Removing frags with length > ", maxFragLength,
        " bp ..\n")
    GA <- GA[width(GA) <= maxFragLength]
    if (length(GA) == 0)
      stop(
        "Fragment filtering resulting in 0 aligned fragments. Please check / change the provided filter size ..\n"
      )
  }
  barcodes <- as.character(GA$barcodeID)
  denom <- table(barcodes)
  uniqueBarcodes <- names(denom)
  id <- factor(barcodes, levels = uniqueBarcodes)
  cat("Finding overlap between peaks and fragments in data ..\n")
  ovPEAKStarts <- findOverlaps(query = peaks,
                               subject = resize(GA,
                                                width = 1, fix = "start"))
  ovPEAKEnds <- findOverlaps(query = peaks,
                             subject = resize(GA,
                                              width = 1, fix = "end"))
  cat(
    "Filtering for valid fragment-peak overlaps based on cut site start/end coordinates ..\n"
  )
  #validHits <- unique.data.frame(rbind(as.data.frame(ovPEAKStarts),as.data.frame(ovPEAKEnds)))
  
  validHits = merge(as.data.table(ovPEAKStarts), as.data.table(ovPEAKEnds),all=T)
  
  cat("Generating matrix of counts ..\n")
  # countdf <-
  #   data.frame(peaks = validHits$queryHits, sample = as.numeric(id)[validHits$subjectHits]) %>%
  #   dplyr::group_by(peaks, sample) %>% dplyr::summarise(count = n()) %>%
  #   data.matrix()
  
  ft = as_fst(data.table(peaks = validHits$queryHits, sample = as.numeric(id)[validHits$subjectHits])) #what happens here? data.table removed from memory?
  
  countdf = ft %>% select_fst(peaks, sample) %>% arrange(peaks, sample) %>% summarise(count = .N, by = list(peaks, sample)) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(
    i = c(countdf[, 1], length(peaks)),
    j = c(countdf[, 2], length(uniqueBarcodes)),
    x = c(countdf[,
                  3], 0)
  )
  colnames(m) <- uniqueBarcodes
  if (addColData) {
    cat("Computing sample read depth and FRIP ..\n")
    colData <-
      data.table(
        sample = uniqueBarcodes,
        depth = as.numeric(denom),
        FRIP = Matrix::colSums(m) / as.numeric(denom),
        stringsAsFactors = FALSE
      )
    stopifnot(all.equal(colData$sample, colnames(m)))
    if (any(colData$FRIP > 1))
      warning(
        "One or more barcodes ended up with FRIP score > 1 .. check your fragment file as it may contain some abnormally large fragments that should be removed ..\n"
      )
    cat("Generating SummarizedExperiment object ..\n")
    SE <-
      SummarizedExperiment(
        rowRanges = peaks,
        assays = list(counts = m),
        colData = colData
      )
  }
  else {
    cat("Generating SummarizedExperiment object ..\n")
    SE <-
      SummarizedExperiment(rowRanges = peaks, assays = list(counts = m))
  }
  cat("Done!\n")
  end_time <- Sys.time()
  cat("Time elapsed: ",
      end_time - start_time,
      units(end_time -
              start_time),
      " \n\n")
  return(SE)
}

getDORCScores = function (ATAC.se, dorcTab, normalizeATACmat = TRUE, geneList = NULL, 
                          nCores = 4) 
{
  if (!all(c("Peak", "Gene") %in% colnames(dorcTab))) 
    stop("The provided gene-peak table must have columns named Peak and Gene ..")
  if (any(dorcTab$Peak > nrow(ATAC.se))) 
    stop("One or more peak indices in the gene-peak table are larger than the total number of peaks in the provided ATAC SE object ..\n Make sure")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% as.character(dorcTab$Gene)))) 
      stop("One or more of the gene names supplied is not present in the gene-peak table provided..\n")
    cat("Running DORC scoring for genes:", geneList, sep = "\n")
    cat("........\n")
    dorcTab <- dorcTab[dorcTab$Gene %in% geneList, ]
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
  }
  else {
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
    cat("Running DORC scoring for all genes in annotation! (n = ", 
        length(dorcGenes), ")\n", sep = "")
  }
  if (normalizeATACmat) {
    cat("Normalizing scATAC counts ..\n")
    ATAC.mat <- assay(centerCounts(ATAC.se, 
                                   chunkSize = 5000))
    gc()
  }
  else {
    cat("Assuming provided scATAC counts are normalized ..\n")
    ATAC.mat <- assay(ATAC.se)
  }
  time_elapsed <- Sys.time()
  cat("Computing DORC scores ..\n")
  cat("Running in parallel using ", nCores, "cores ..\n")
  dorcMatL <- pbmcapply::pbmclapply(X = dorcGenes, FUN = function(x) {
    dorcPeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% x])
    if (length(dorcPeaks) > 1) {
      dorcCounts <- Matrix::colSums(ATAC.mat[dorcPeaks, 
      ])
    }
    else if (length(dorcPeaks == 1)) {
      dorcCounts <- ATAC.mat[dorcPeaks, ]
    }
  }, mc.cores = nCores)
  dorcMat <- Matrix(do.call("rbind", dorcMatL), sparse = TRUE)
  rownames(dorcMat) <- dorcGenes
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)), 
      "\n\n")
  return(dorcMat)
}