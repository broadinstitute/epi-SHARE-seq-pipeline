output <- ''
genome <- '' # hg19, mm10, etc
bcFile <- '' # barcodes.txt
readFile <- '' #filtered.counts.csv.gz
h5file <- '' #h5 UMI counts
qc_file <- '' #counts.csv.gz

ncores <- 4
read_thresh <- 300
umi_thresh <- 500
gtfPATH <- '/mnt/users/snagaraja/ref/gtf/'

### VERSION CONTROL
## v2.1: reformatted to pull functions from separate file
## v2.0: adjusted for SHARE-seq V4.12 pipeline output

myPath <- .libPaths()
myPath <- c(myPath,'/mnt/bin/R/library/')
.libPaths(myPath)
source('/home/snagaraja/comp_tools/scRNA_scATAC_scripts/function_code/barcode_filtering.v1.1.R')
library(pbmcapply)
library(hdf5r)

bc <- read.csv(bcFile,header=F)
readCounts <- read.csv(readFile)
qc <- read.csv(qc_file,header=T)
umiCounts <- ReadCB_h5(h5file, use.names = F, unique.features = F)

getGTF <- getGtfGenes(gtfPATH,genome)

bc$bcID <- with(bc,paste(V1,V2,V3,V4,sep=','))
readCounts$bcID <- with(readCounts,paste(R1,R2,R3,P5,sep=','))
readCounts <- readCounts[(readCounts$bcID %in% bc$bcID),]

## filter by total fragments per barcode
pdf(paste0(output,'.frag_distribution.pdf'),width=5,height=4)
with(readCounts,hist(log10(fragments+1),main='Fragments per barcode'))
abline(v=log10(read_thresh+1),lty=2,col='red')
dev.off()
bc_pf1 <- readCounts$bcID[(readCounts$fragments >= read_thresh)]

colnames(qc)[grepl('_libsize',colnames(qc))] <- 'libsize'
pdf(paste0(output,'.lib_size_pf1.pdf'),width=5,height=5)
qc_pf1 <- qc[(with(qc,paste(R1,R2,R3,P5,sep=',')) %in% bc_pf1),]
qc_pf1 <- qc_pf1[order(qc_pf1$libsize,decreasing=T),]
plot(log10(qc_pf1$libsize),ylab='log10(Library size)',xlab='Barcode rank',main='Estimated Lib Size - PF1',pch=16,col='blue')
avg <- floor(median(qc_pf1$libsize))
text(0.9*nrow(qc_pf1),0.9*max(log10(qc_pf1$libsize)),labels=paste0('Median lib size: ',avg),pos=2)
dev.off()

## check barcode wells
wells <- wellReadPlot(bc_pf1,readCounts,3)
pdf(paste0(output,'.frag_well_distribution.pf1.pdf'),width=5,height=3)
plotShareWells(wells[[1]],1)
plotShareWells(wells[[2]],2)
plotShareWells(wells[[3]],3)
dev.off()

## filter by genes detected by barcodes
umi_pf1 <- as(as.matrix(umiCounts[,(colnames(umiCounts) %in% bc_pf1)]),'sparseMatrix')
totUMI <- data.frame(bc=colnames(umi_pf1),totUMI=Matrix::colSums(umi_pf1))

pdf(paste0(output,'.pf1_total_UMIs.pdf'),width=5,height=4)
hist(log10(totUMI$totUMI+1),main='UMIs per barcode',xlab='log10(Total UMIs)',breaks=40)
abline(v=log10(umi_thresh+1),lty=2,col='red')
dev.off()

umi_pf2 <- umi_pf1[,(totUMI$totUMI > umi_thresh)]
if (getGTF[[1]]){
  cat('Adding in missing genes from GTF to UMI matrix\n')
  genes <- getGTF[[2]]
  missing_genes <- genes[!(genes %in% rownames(umi_pf2))]
  rest_mat <- matrix(nrow=length(missing_genes),ncol=ncol(umi_pf2),data=0)
  rownames(rest_mat) <- missing_genes
  combined <- rbind(as.matrix(umi_pf2),rest_mat)
  combined <- combined[match(genes,rownames(combined)),]
  umi_pf2 <- as(combined,'sparseMatrix')
}
saveRDS(umi_pf2,paste0(output,'.UMIcounts.pf2.rds'))
writeLines(colnames(umi_pf2),paste0(output,'.bc_pf2.txt'))

#check wells again
bc_pf2 <- colnames(umi_pf2)

pdf(paste0(output,'.lib_size_pf2.pdf'),width=5,height=5)
qc_pf2 <- qc[(with(qc,paste(R1,R2,R3,P5,sep=',')) %in% bc_pf2),]
qc_pf2 <- qc_pf2[order(qc_pf2$libsize,decreasing=T),]
plot(log10(qc_pf2$libsize),ylab='log10(Library size)',xlab='Barcode rank',main='Estimated Lib Size - PF2',pch=16,col='blue')
avg <- floor(median(qc_pf2$libsize))
text(0.9*nrow(qc_pf2),0.9*max(log10(qc_pf2$libsize)),labels=paste0('Median lib size: ',avg),pos=2)
dev.off()

wells2 <- wellReadPlot(bc_pf2,readCounts,3)
pdf(paste0(output,'.bc_well_distribution.pf2.pdf'),width=5,height=3)
plotShareWells(wells2[[1]],1)
plotShareWells(wells2[[2]],2)
plotShareWells(wells2[[3]],3)
dev.off()
