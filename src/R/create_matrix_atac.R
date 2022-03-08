#!/usr/bin/Rscript
# generate atac count matrix

args <- commandArgs(); 
print(args)

suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(DropletUtils))
suppressMessages(library(GenomicRanges))
suppressMessages(library(SummarizedExperiment))

fragFile = args[6] #ATAC fragment file (bed, bedpe)
peakFile = args[7] #Peak file (cCRE)
outputFile = args[8] #Output filename, h5 file with genomic peaks (rows) and cell barcodes (columns)

getCountsFromFrags <- function (fragRanges,
                                peaks,
                                barcodeList = NULL,
                                maxFragLength = NULL,
                                addColData = TRUE) {
  start_time <- Sys.time()
  
  GA <- fragRanges
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
  validHits <- unique.data.frame(rbind(as.data.frame(ovPEAKStarts),
                                       as.data.frame(ovPEAKEnds)))
  require(dplyr)
  cat("Generating matrix of counts ..\n")
  countdf <-
    data.frame(peaks = validHits$queryHits, sample = as.numeric(id)[validHits$subjectHits]) %>%
    dplyr::group_by(peaks, sample) %>% dplyr::summarise(count = n()) %>%
    data.matrix()
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
      data.frame(
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


#peaks <- read.table("data/GM_nonoverlap.bed",sep="\t", header=F)
#frags <- read.table("data/shareseq-project.atac.GRCh38.cleaned.filtered.bedpe",sep="\t", header=F)

peaks <- read.table(peakFile,sep="\t", header=F)
frags <- read.table(fragFile,sep="\t", header=F)

peakRanges <- makeGRangesFromDataFrame(peaks,seqnames.field = "V1",start.field = "V2",end.field = "V3",starts.in.df.are.0based = TRUE)

fragRanges = makeGRangesFromDataFrame(frags, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

# Get counts from atac frag file
peaksSE <- getCountsFromFrags(fragRanges=fragRanges, peaks=peakRanges)

head(peaksSE)
names <- paste0(peaks$V1, ':', peaks$V2, '-', peaks$V3)

write10xCounts(
  path = outputFile,
  x = assay(peaksSE),
  barcodes = colnames(peaksSE),
  gene.id = names,
  gene.symbol = names,
  overwrite = FALSE
)