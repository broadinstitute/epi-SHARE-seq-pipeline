output <- '' # output prefix
genome <- ''
fragBed <- '' #raw fragments.tsv.gz file
bcFile <- '' #bc_pf from refilter_barcodes
peakFile <- '' #fixedwidth and resized peaks
qcFile <- '' #counts.csv.gz
# both declump and doublet cell barcodes to exclude
rna_out_files <- c('',
                   ''
)

#for plotting genes on UMAP
rna_smooth_rds <- ''
markerList <- readLines('~/comp_tools/scRNA_scATAC_scripts/gene_lists/')
markerList <- append(markerList,c(''))

## for TF - RNA cor plot
tfToAdd <- c('')

# for motif score plotting
TFsToPlot <- c('')

## Ad barcode matching, c('P1.RNA','P1.ATAC')
isMultiSample <- T #T if more than one Ad1.xx per RNA/ATAC
adList <- list( c('',''),
                c('','')
)


### Version Control
# v2.6: added by passage plotting feature, made gene scores boolean
# v2.1: added RNA-ATAC Ad1 matching for large datasets where collisions outside Ad1.xx may occur; pull functions by sourcing files
# v2.0: adapted for outputs of SHARE-seq V4.12 pipeline
# v1.3: updated cisTopic with R upgrade on workstation

size.thresh <- 500
frip.thresh <- 0.2
num_topics <- 80
umap_pt_size <- 0.5
zCap <- 3
cores <- 4
nBg <- 250
nTFplot <- 30
nTop <- 50
doGeneScores <- F

###

myPath <- .libPaths()
myPath <- c(myPath,'/mnt/bin/R/library/')
.libPaths(myPath)
source('/home/snagaraja/comp_tools/scRNA_scATAC_scripts/function_code/singlet_clustering.shared.v1.6.R')
source('/home/snagaraja/comp_tools/scRNA_scATAC_scripts/function_code/singlet_clustering.ATAC.v1.6.R')

### General QC, filtering 

## Getting reads in peaks
peaks <- read.table(peakFile,stringsAsFactors=F)
colnames(peaks) <- c('chr','start','end')
peaks <- makeGRangesFromDataFrame(peaks)
bcList <- readLines(bcFile)
if (ncol(stri_split_fixed(bcList,'.',simplify=T)) == 8){
  cell_split <- data.frame(stri_split_fixed(bcList,'.',simplify=T))
  bcList <- with(cell_split,paste0(X1,'.',X2,',',X3,'.',X4,',',X5,'.',X6,',',X7,'.',X8))
}

fragments <- Signac::CreateFragmentObject(path = fragBed, cells = bcList, validate.fragments = FALSE)
fm <- Signac::FeatureMatrix(
  fragments,
  peaks,
  cells = bcList,
  verbose = TRUE
)
se <- SummarizedExperiment(list(counts=fm))
rowRanges(se) <- peaks

## Filtering cells by QC
qc <- read.csv(qcFile)

qc$barcode <- with(qc,paste(R1,R2,R3,P5,sep=','))
qc_sub <- qc[(qc$barcode %in% colnames(se)),]
colnames(qc_sub)[grepl('_libsize',colnames(qc_sub))] <- 'libsize'
colData(se)$libsize <- qc_sub$libsize[match(colnames(se),qc_sub$barcode)]
colData(se)$tot_reads <- qc_sub[,grepl('_unique',colnames(qc_sub))]
colData(se)$depth <- colSums(assays(se)$counts)
colData(se)$FRIP <- with(colData(se),(depth/tot_reads))

saveRDS(se,paste(output,'all_counts.rds',sep='.'))

pf <- plotAtacQC(se,var1='libsize',thresh1=size.thresh,frip.thresh=frip.thresh)
cat('Library size cutoff:',size.thresh,'\n')
cat('FRIP cutoff:',frip.thresh,'\n')
cont <- readline(prompt='Change cutoffs? (y/n): ')

while (cont=='y' | cont=='Y'){
  size.thresh <- as.numeric(readline(prompt='New size threshold: '))
  frip.thresh <- as.numeric(readline(prompt='New FRIP threshold: '))
  pf <- plotAtacQC(se,var1='libsize',thresh1=size.thresh,frip.thresh=frip.thresh)
  cont <- readline(prompt='Change cutoffs? (y/n): ')
}
dev.off()
pdf(paste0(output,'.qc_plot.pdf'),width=4,height=4)
pf <- plotAtacQC(se,var1='libsize',thresh1=size.thresh,frip.thresh=frip.thresh)
dev.off()

se.qc <- se[,pf]
saveRDS(se.qc,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.all.rds'))
rm(se,fm)
gc()

#############################################################
################# 01 - ATAC ONLY DECLUMPING #################
#############################################################

## all cell cisTopics
cis.all <- getCisTopics(se.qc,'all',num_topics)
topics <- cis.all[[1]]
write.table(topics,paste0(output,'.cisTopics_topics_zCell.all.txt'),sep='\t',quote=F)
write.table(cis.all[[2]],paste0(output,'.cisTopics_region_data.all.txt'),sep=',',quote=F)

## UMAP
UMAP <- cacheUMAP(	cachePath = './umapCache',
                   dataSet = topics,
                   seed = 123)
umap.t <- data.frame(UMAP$layout)
colnames(umap.t) <- c('UMAP1','UMAP2')
saveRDS(umap.t,paste0(output,'.cisTopic_',ncol(topics),'.umap_embed.all.rds'))
rm(UMAP)
gc()

## Finding KNN
knn_data <- get.knn(topics,k=20)
knn <- knn_data$nn.index

# QC filtering
qc_sub.all <- qc_sub[(qc_sub$barcode %in% rownames(topics)),]
qc_sub.all <- cbind(qc_sub.all,umap.t)

pdf(paste0(output,'.umap.log10libsize.all.pdf'),width=4.5,height=4)
g <- ggplot(qc_sub.all,aes(x=UMAP1,y=UMAP2,color=log10(libsize))) + geom_point(size=1,stroke=0) + 
  scale_color_gradient(low='navy',high='yellow') +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()

clump.ind <- getShareClumps(qc_sub.all$libsize,knn,umap.t,ptSize=1)
writeLines(rownames(umap.t)[clump.ind],paste0(output,'.ATAC_clump_barcodes.txt'))
se.d <- se.qc[,!(clump.ind)]
saveRDS(se.d,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.declump.rds'))

#############################################################
################ 02 - FILTERING RNA BARCODES ################
#############################################################

rna_bc_out_list <- vector("list",length=length(rna_out_files))
for (i in 1:length(rna_out_files)){
  rna_bc_out_list[[i]] <- readLines(rna_out_files[i])
}
rna_bc_out <- unique(unlist(rna_bc_out_list))

matchBC <- matchShareSamples( atac_bc = colnames(se.d),
                              rna_bc = rna_bc_out,
                              isMulti = isMultiSample,
                              adList = adList
                              )
keep_ind <- !(colnames(se.d) %in% matchBC[[1]])
se.s <- se.d[,keep_ind]
saveRDS(se.s,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.singlet.rds'))

#############################################################
################### 03 - SINGLET ANALYSIS ###################
#############################################################

cis.s <- getCisTopics(se.s,'singlet',nTopics = seq(ncol(topics)-10,ncol(topics)+10,by=10))
topics.s <- cis.s[[1]]
write.table(topics.s,paste0(output,'.cisTopics_topics_zCell.singlet.txt'),sep='\t',quote=F)
write.table(cis.s[[2]],paste0(output,'.cisTopics_region_data.singlet.txt'),sep=',',quote=F)

## UMAP
UMAP.s <- cacheUMAP(	cachePath = './umapCache',
                   dataSet = topics.s,
                   seed = 123)
umap.s <- data.frame(UMAP.s$layout)
colnames(umap.s) <- c('UMAP1','UMAP2')
saveRDS(umap.s,paste0(output,'.cisTopic_',ncol(topics.s),'.umap_embed.singlet.rds'))

## Finding KNN
knn_data.s <- get.knn(topics.s,k=20)
knn.s <- knn_data.s$nn.index
write.table(knn.s,paste0(output,'.cisTopics_knn_20.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')

# plot depth
qc_sub.s <- qc_sub[(qc_sub$barcode %in% rownames(topics.s)),]
qc_sub.s$FRIP <- colData(se.s)$FRIP

## if multiple samples per run
#r1.s <- as.numeric(stri_split_fixed(qc_sub.s$R1,'.',simplify=T)[,2])
#qc_sub.s$sample <- 'blank'
#qc_sub.s$rep[(r1.s<=24)] <- 'pbs6male'
r1.s <- as.numeric(stri_split_fixed(qc_sub.s$R1,'.',simplify=T)[,2])
qc_sub.s$sample <- 'blank'
qc_sub.s$replicate <- 'blank'
qc_sub.s$batch <- 'blank'

batch_p5 <- (qc_sub.s$P5 %in% c('P1.64','P1.60','P1.61'))
qc_sub.s$batch[batch_p5] <- 'DSS04'
qc_sub.s$sample[batch_p5 &(r1.s<=32)] <- 'DSS_P1'
qc_sub.s$sample[batch_p5 &(r1.s>32) & (r1.s<=48)] <- 'Control_P1'
qc_sub.s$sample[batch_p5 &(r1.s>48) & (r1.s<=80)] <- 'DSS_P4'
qc_sub.s$sample[batch_p5 &(r1.s>80) & (r1.s<=96)] <- 'Control_P4'

qc_sub.s$replicate[batch_p5 & (r1.s<=12)] <- 'DSS04_m1_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>12) & (r1.s<=24)] <- 'DSS04_m2_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>24) & (r1.s<=32)] <- 'DSS04_m3_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>32) & (r1.s<=40)] <- 'DSS04_m4_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>40) & (r1.s<=48)] <- 'DSS04_m5_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>48) & (r1.s<=60)] <- 'DSS04_m1_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>60) & (r1.s<=72)] <- 'DSS04_m2_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>72) & (r1.s<=80)] <- 'DSS04_m3_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>80) & (r1.s<=88)] <- 'DSS04_m4_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>88) & (r1.s<=96)] <- 'DSS04_m5_P4'


batch_p5 <- (qc_sub.s$P5 %in% c('P1.53','P1.62','P1.57'))
qc_sub.s$batch[batch_p5] <- 'DSS07'
qc_sub.s$sample[batch_p5 & (r1.s<=24)] <- 'DSS_P1'
qc_sub.s$sample[batch_p5 & (r1.s>24) & (r1.s<=48)] <- 'Control_P1'
qc_sub.s$sample[batch_p5 & (r1.s>48) & (r1.s<=72)] <- 'DSS_P4'
qc_sub.s$sample[batch_p5 & (r1.s>72) & (r1.s<=96)] <- 'Control_P4'

qc_sub.s$replicate[batch_p5 & (r1.s<=12)] <- 'DSS07_m1_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>12) & (r1.s<=24)] <- 'DSS07_m2_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>24) & (r1.s<=36)] <- 'DSS07_m3_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>36) & (r1.s<=48)] <- 'DSS07_m4_P1'
qc_sub.s$replicate[batch_p5 & (r1.s>48) & (r1.s<=60)] <- 'DSS07_m1_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>60) & (r1.s<=72)] <- 'DSS07_m2_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>72) & (r1.s<=84)] <- 'DSS07_m3_P4'
qc_sub.s$replicate[batch_p5 & (r1.s>84) & (r1.s<=96)] <- 'DSS07_m4_P4'

write.table(qc_sub.s,paste0(output,'.qc_stats.singlet.txt'),quote=F,sep='\t',row.names=F)

qc_sub.s <- cbind(qc_sub.s,umap.s)
pdf(paste0(output,'.umap.log10libsize.singlet.pdf'),width=4.5,height=4)
g <- ggplot(qc_sub.s,aes(x=UMAP1,y=UMAP2,color=log10(libsize))) + geom_point(size=1,stroke=0) + 
  scale_color_gradient(low='navy',high='yellow') +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()

rm(UMAP.s,knn_data.s,g)
gc()

if (sum(colnames(qc_sub.s) == 'sample') > 0){
  pdf(paste0(output,'.umap_singlet.sample_knn_density.pdf'),width=5,height=4)
  for (S in unique(qc_sub.s$sample)){
    cat(S,'\n')
    plotGroupDensity( G = S,
                      U = umap.s,
                      groups = qc_sub.s$sample,
                      NN = knn.s,
                      ncores = cores
    )
  }
  dev.off()
}

if (sum(colnames(qc_sub.s) == 'phase') > 0){
  pdf(paste0(output,'.umap_singlet.phase_knn_density.pdf'),width=5,height=4)
  for (P in unique(qc_sub.s$phase)){
    cat(P,'\n')
    plotGroupDensity( G = P,
                      U = umap.s,
                      groups = qc_sub.s$phase,
                      NN = knn.s,
                      ncores = cores
    )
  }
  dev.off()
}

if (sum(colnames(qc_sub.s) == 'batch') > 0){
  pdf(paste0(output,'.umap_singlet.batch_knn_density.pdf'),width=5,height=4)
  for (B in unique(qc_sub.s$batch)){
    cat(B,'\n')
    plotGroupDensity( G = B,
                      U = umap.s,
                      groups = qc_sub.s$batch,
                      NN = knn.s,
                      ncores = cores
    )
  }
  dev.off()
}

#louvain
ig <- graph.empty(nrow(knn.s))
edge_list <- pbmclapply(X=1:nrow(knn.s),FUN=function(x){
  return(unlist(mapply(c,rep(x,ncol(knn.s)),knn.s[x,],SIMPLIFY=F)))
},mc.cores=cores)
ig <- add_edges(ig,unlist(edge_list))
comm <- cluster_louvain(as.undirected(ig))
clust <- data.frame(sample=rownames(topics.s),community=as.vector(membership(comm)),stringsAsFactors = F)
write.table(clust,paste(output,'cisTopics_singlet.louvain_clusters.txt',sep='.'),quote=F,sep='\t',row.names=F)

u <- cbind(umap.s,as.character(clust$community))
colnames(u)[3] <- 'community'

pdf(paste0(output,'.cisTopics_umap_singlet.louvain_clusters.pdf'),width=5,height=4)
g <- ggplot(u,aes(x=UMAP1,y=UMAP2,color=community)) + geom_point(size=umap_pt_size,stroke=0) + 
  scale_color_viridis(discrete = T, option = "D")+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
for (comm in unique(u$community)){
  g <- ggplot(u,aes(x=UMAP1,y=UMAP2,color=(community==comm))) + geom_point(size=umap_pt_size,stroke=0) +
    guides(color=guide_legend(title=paste0('Community ',comm))) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g)
}
dev.off()

rm(edge_list,g,ig)
gc()

if (sum(qc_sub.s$sample != 'blank') > 0){
  plotBySampleATAC(qc_sub.s,clust)
}

if (doGeneScores){
  # gene scores 
  gs <- getGeneScoresFromFrags( fragFile = fragBed,
                                barcodeList = colnames(se.s),
                                genome=genome,
                                decayMaxFrac = 0.01)
  saveRDS(gs,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.gene_scores_frags.singlet.rds'))
  
  gs.smooth <- smoothCellsNN( E = gs,
                              K = knn.s,
                              cores = cores
  )
  saveRDS(gs.smooth,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.gene_scores_frags.smooth.singlet.rds'))
}

# plot genes on UMAP
rna.smooth <- readRDS(rna_smooth_rds)
matchBC.s <- matchShareSamples( atac_bc = colnames(se.s),
                              rna_bc = colnames(rna.smooth),
                              isMulti = isMultiSample,
                              adList = adList
                              )
cat('Total ATAC barcodes:',ncol(se.s),'\n')
cat('Total RNA barcodes:',ncol(rna.smooth),'\n')
cat('Shared barcodes between ATAC and RNA:',length(matchBC.s[[1]]),'\n')

rna.sub <- rna.smooth[,matchBC.s[[2]]]
umap.sub <- umap.s[matchBC.s[[1]],]

if (!dir.exists('./umap_gene_plots')){
  dir.create('./umap_gene_plots')
}
setwd('./umap_gene_plots')
for (marker in markerList){
  if (doGeneScores){
    pdf(paste0(output,'.umap_singlet.gene_score.',marker,'.pdf'),width=4.5,height=4)
    plotUMIgene(umap.s,gs.smooth,marker,zCap,umap_pt_size)
    dev.off()
  }
  
  pdf(paste0(output,'.umap_singlet.RNA.',marker,'.pdf'),width=4.5,height=4)
  plotUMIgene(umap.sub,rna.sub,marker,zCap,umap_pt_size)
  dev.off()
}
setwd('../')

# motif scores
cat('Now computing motif scores with chromVAR\n')
BiocParallel::register(BiocParallel::MulticoreParam(cores, progressbar = TRUE))
if (genome == 'hg19'){
  library(BSgenome.Hsapiens.UCSC.hg19)
  BSg <-BSgenome.Hsapiens.UCSC.hg19
  motifSet <- human_pwms_v2
  organism <- 'human'
} else if (genome == 'mm10'){
  library(BSgenome.Mmusculus.UCSC.mm10)
  BSg <- BSgenome.Mmusculus.UCSC.mm10
  motifSet <- mouse_pwms_v3
  organism <- 'mouse'
} else if (genome == 'hg38'){
  library(BSgenome.Hsapiens.UCSC.hg38)
  BSg <- BSgenome.Hsapiens.UCSC.hg38
  motifSet <- human_pwms_v2
  organism <- 'human'
} else {
  cat('Genome not recognized')
}
se.s <- addGCBias(se.s,genome=BSg)
se.filt <- filterPeaks(se.s)
saveRDS(se.filt,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.singlet_peakFilter.rds'))
rm(se.s,gs,gs.smooth,rna.smooth,rna.sub)
gc()

set.seed(123)
bg <- getBackgroundPeaks(se.filt,niterations=nBg)
motif_ix <- matchMotifs(pwms = motifSet,subject = se.filt,genome=BSg)
dev_motif <- computeDeviations(object = se.filt,annotations = motif_ix,background_peaks=bg)
saveRDS(dev_motif,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.dev_motif.rds'))
rm(bg)
gc()

## plot RNA - motif correlation
rownames(dev_motif) <- extractTFNames(rownames(dev_motif))
z <- assays(dev_motif)$z

pdf(paste0(output,'.TFmotif_expr_cor.pdf'),width=5,height=5)
d <- tfExprCor( scores = z,
                expr = rna.smooth,
                multiSample = isMultiSample,
                adList = adList,
                addGenes = tfToAdd,
                nLabel = nTFplot,
                cores = cores
                )
dev.off()
write.table(d,paste0(output,'.TFmotif_expr_cor.csv'),sep=',',row.names=F,quote=F)

# plot a few motif scores on UMAP
TFsToPlot <- unique(append(TFsToPlot,c('Tead4')))
pdf(paste0(output,'.umap_singlet.TF_scores.pdf'),width=4.5,height=4)
for (tf in TFsToPlot){
  cat(tf,'\n')
  plotUMIgene(umap.s,z,tf,zCap,umap_pt_size,palette='solar_extra')
}
dev.off()

### by passage
dir.create('./by_passage/')
setwd('./by_passage/')
passages <- c('P1','P4')
for (passage in passages){
  cat(passage,'\n')
  passage_ind <- grepl(passage,qc_sub.s$sample)
  qc_p <- qc_sub.s[passage_ind,]
  topics_scores_p <- topics.s[passage_ind,]
  
  UMAP.pass <- cacheUMAP(	cachePath = './umapCache',
                          dataSet = topics_scores_p,
                          seed = 123)
  umap.pass <- data.frame(UMAP.pass$layout)
  colnames(umap.pass) <- c('UMAP1','UMAP2')
  saveRDS(umap.pass,paste0(output,'.passage_',passage,'.umap_coord.singlet.rds'))
  rm(UMAP.pass)
  gc()
  
  # Smoothening KNN
  knn_data.pass <- get.knn(topics_scores_p,k=20)
  knn.pass <- knn_data.pass$nn.index
  write.table(knn.pass,paste0(output,'.passage_',passage,'.knn_20.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')
  
  ig <- graph.empty(nrow(knn.pass))
  edge_list <- pbmclapply(X=1:nrow(knn.pass),FUN=function(x){
    return(unlist(mapply(c,rep(x,ncol(knn.pass)),knn.pass[x,],SIMPLIFY=F)))
  },mc.cores=cores)
  ig <- add_edges(ig,unlist(edge_list))
  comm <- cluster_louvain(as.undirected(ig))
  clust.pass <- data.frame(sample=rownames(umap.pass),community=as.vector(membership(comm)),stringsAsFactors = F)
  write.table(clust.pass,paste0(output,'.passage_',passage,'.cisTopics_singlet.louvain_clusters.txt'),quote=F,sep='\t',row.names=F)
  
  u.pass <- cbind(umap.pass,as.character(clust.pass$community))
  colnames(u.pass)[3] <- 'community'
  
  pdf(paste0(output,'.passage_',passage,'.cisTopics_umap_singlet.louvain_clusters.pdf'),width=5,height=4)
  g <- ggplot(u.pass,aes(x=UMAP1,y=UMAP2,color=community)) + geom_point(size=umap_pt_size,stroke=0) + 
    scale_color_viridis(discrete = T, option = "D")+
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g)
  for (comm in unique(u.pass$community)){
    g <- ggplot(u.pass,aes(x=UMAP1,y=UMAP2,color=(community==comm))) + geom_point(size=umap_pt_size,stroke=0) +
      guides(color=guide_legend(title=paste0('Community ',comm))) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
    print(g)
  }
  dev.off()
  
  rm(edge_list,g,ig)
  gc()
  
  qc_p$UMAP1 <- umap.pass$UMAP1
  qc_p$UMAP2 <- umap.pass$UMAP2
  
  if (sum(colnames(qc_p) == 'sample') > 0){
    pdf(paste0(output,'.umap_singlet.passage_',passage,'.sample_knn_density.pdf'),width=5,height=4)
    for (S in unique(qc_p$sample)){
      cat(S,'\n')
      plotGroupDensity( G = S,
                        U = umap.pass,
                        groups = qc_p$sample,
                        NN = knn.pass,
                        ncores = cores
      )
    }
    dev.off()
  }
  if (sum(colnames(qc_p) == 'batch') > 0) {
    pdf(paste0(output,'.umap_singlet.passage_',passage,'.by_batch.pdf'),width=5,height=4)
    g1 <- ggplot(shuf(qc_p),aes(x=UMAP1,y=UMAP2,color=batch)) + geom_point(size=umap_pt_size,stroke=0) +
      scale_color_viridis(discrete = T, option = "D") +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
    print(g1)
    dev.off()
  }
  
  #plot some genes
  markerList2 <- c('Lgr5','Clu','Ly6a')
  umap.pass.sub <- umap.pass[(rownames(umap.pass) %in% rownames(umap.sub)),]
  rna.pass.sub <- rna.sub[,(rownames(umap.sub) %in% rownames(umap.pass.sub))]
  pdf(paste0(output,'.passage_',passage,'.umap_singlet.markers.pdf'),width=4.5,height=4)
  for (marker in markerList2){
    cat(marker,'\n')
    plotUMIgene(umap.pass.sub,rna.pass.sub,marker,zCap,2*umap_pt_size)
  }
  dev.off()
  
  #plot some motifs
  z.pass <- z[,passage_ind]
  pdf(paste0(output,'.passage_',passage,'.umap_singlet.TF_scores.pdf'),width=4.5,height=4)
  for (tf in TFsToPlot){
    cat(tf,'\n')
    plotUMIgene(umap.pass,z.pass,tf,zCap,umap_pt_size,palette='solar_extra')
  }
  dev.off()
}


cat('Done!\n')
