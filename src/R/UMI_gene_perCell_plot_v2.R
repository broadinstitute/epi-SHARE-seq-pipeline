#!/usr/bin/Rscript
# plot umi/cell or gene/cell

args <- commandArgs(); # print(args)

suppressMessages(library(matrixStats))
suppressMessages(library(reshape2))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressPackageStartupMessages(library("DropletUtils"))
suppressPackageStartupMessages(library("rhdf5"))
suppressPackageStartupMessages(library("SummarizedExperiment"))

# load and convert to count matrix
outFile1 <- "out.UMIcounts.csv.gz"
#cutoff.bed.gz
inFile1 <- args[6]

print("Loading linear gene table...")
linear <- fread(inFile1, header = F, sep="\t")
print("Finished loading file")

# Df2 <- acast(linear[1:10000, ], V1~V2, value.var="V3", fill=0)
# Df2 <- acast(linear, V2~V1, value.var="V3", fill=0)
# Df2 <- as(acast(linear, V2~V1, value.var="V3", fill=0), "sparseMatrix")
# head(temp)

cells <- unique(linear$V1)
genes <- unique(linear$V2)
step=20000

print(paste("No. genes x cells = ", length(genes), "x", length(cells), sep=""))

if (length(cells) > step){
   window <- data.frame(start=seq(1, length(cells), step), end=c(seq(step, length(cells), step), length(cells)))
   for (i in 1:nrow(window)){
      print(window[i, ])
      sub <- linear %>% filter(V1 %in% cells[window[i,1]:window[i,2]])
      Mx.sub <- as(acast(sub, V2~V1, value.var="V3", fill=0), "sparseMatrix")
      emptygenes <- genes[!genes %in% rownames(Mx.sub)]
      emptyMx <- as(matrix(data = 0, ncol = ncol(Mx.sub), nrow = length(emptygenes)), "sparseMatrix")
      colnames(emptyMx) <- colnames(Mx.sub)
      rownames(emptyMx) <- emptygenes
      Mx.sub <- rbind(Mx.sub,emptyMx)
      if (i ==1){
        Df2 <- Mx.sub
      } else{
        Df2 <- cbind(Df2, Mx.sub[rownames(Df2), ])
      }
   }
} else {
  Df2 <- as(acast(linear, V2~V1, value.var="V3", fill=0), "sparseMatrix")
}

# print(paste("Output to ", outFile1, sep=""))
# write.table(as.matrix(Df2), gzfile(outFile1), sep = "\t", col.names=T, row.names=T, quote=F)

## make h5
suppressPackageStartupMessages(library("DropletUtils"))
suppressPackageStartupMessages(library("rhdf5"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("Matrix"))
rna.count <- Df2
print(paste("Writing to ", "out.gene.bc.matrices.h5", sep=""))


write10xCounts(
  path = "out.gene.bc.matrices.h5",
  x = rna.count,
  barcodes = colnames(rna.count),
  gene.id = rownames(rna.count),
  gene.symbol = rownames(rna.count),
  overwrite = FALSE
)

print("Finished making h5 file for sample")

## make plots
if (ncol(Df2) <= 50000){
   Df2 <- as.matrix(Df2)
   umiSum <- colSums(Df2);
   geneSum <- as.data.frame(colCounts(Df2>0))
   rownames(geneSum) <- colnames(Df2)

   print("top10 UMI counts: ")
   head(sort(umiSum, decreasing = T), 10)

   # set cutoff of gene
   Cutoff <- 0
   Idx <- geneSum > Cutoff; # sum(Idx) # number of cell pass filter
   # mean(geneSum[Idx])
   colnames(geneSum) <- "Genes"
   Df3 <- as.data.frame(sort(geneSum$Genes, decreasing = T))
   Df3$Count <- c(1: nrow(Df3)); colnames(Df3) <- c("Gene", "Count")

   file3 <- 'out.detected.genes.pdf'
   pdf(file3)
plot(Df3$Count,log10(Df3$Gene), xlab="Barcode rank", ylab = "log10 (Genes)", main = "Detected Genes per Cell", col="darkblue", pch=16)
garbage <- dev.off()

# set cutoff of UMI
Cutoff2 <- 0
Idx2 <- umiSum > Cutoff2; # sum(Idx) # number of cell pass filter
# mean(umiSum[Idx2])
Df4 <- as.data.frame(sort(umiSum, decreasing = T))
Df4$Count <- c(1: nrow(Df4)); colnames(Df4) <- c("Umi", "Count")

## UMIs per cell
# head(umiSum); length(umiSum)
umiSum <- as.data.frame(umiSum)
geneSum <- as.data.frame(geneSum)

# plot umi vs genes
Df <- cbind(geneSum$Genes, umiSum[match(rownames(geneSum), rownames(umiSum)), ]); colnames(Df) <- c("Genes","UMIs")
# head(Df)
Df <- as.data.frame(Df)

file7 <- 'out.UMIvsGenes.pdf'
pdf(file7)
plot(Df$UMIs,Df$Genes, xlab="UMIs", ylab = "Detected Genes", main = "Detected Genes per Cell", col="darkblue", pch=16)
legend("topleft", c(paste("Number of Cells that detected > 0 genes:       ",sum(Df$Genes > 0),sep = ""),
paste("Number of Cells that detected > 10 genes:     ",sum(Df$Genes > 10), sep = ""),
paste("Number of Cells that detected > 100 genes:   ",sum(Df$Genes > 100),sep = ""),
paste("Number of Cells that detected > 500 genes:   ",sum(Df$Genes > 500),sep = ""),
paste("Number of Cells that detected > 1000 genes: ",sum(Df$Genes > 1000), sep = "")), bty="n")
garbage <- dev.off()

print(paste("Number of Cells that detected > 0 genes:    ",sum(Df$Genes > 0), sep = ""))
print(paste("Number of Cells that detected > 10 genes:   ",sum(Df$Genes > 10), sep = ""))
print(paste("Number of Cells that detected > 100 genes:  ",sum(Df$Genes > 100), sep = ""))
print(paste("Number of Cells that detected > 500 genes:  ",sum(Df$Genes > 500), sep = ""))
print(paste("Number of Cells that detected > 1000 genes: ",sum(Df$Genes > 1000), sep = ""))
} else {
  print("More than 50k barcodes detected, skip ploting Gene_per_cell...")
}
