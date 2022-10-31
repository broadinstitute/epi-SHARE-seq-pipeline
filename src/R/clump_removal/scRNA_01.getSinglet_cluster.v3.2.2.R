output <- '' #prefix for output files
umiFile <- '' #pf2.rds UMI counts
qc_file <- '' #counts.csv.gz file

two_thresh <- F #use two size thresholds, eg one for lymphocytes and one for epithelial
if (two_thresh){
  thresh_marker <- 'Ptprc'
  umi_thresh1 <- 300
  umi_thresh2 <- 800
}

isMultiSample <- T # if multiple samples per run, specify further in file and will plot UMAP / QC by sample
cores <- 4
umi_thresh <- 1000
doub_est <- 0.06
umap_pt_size <- 0.5
zCap <- 3
nTop <- 50
markerList <- readLines('~/comp_tools/scRNA_scATAC_scripts/gene_lists/')

## Version Control
# v3.2 / 3.2.1: sped up louvain clustering, adjusted smoothening function for large datasets
# v3.1: plotting sample density UMAP
# v3.0: added feature to filter twice on size to allow for cell types with less RNA
# v2.1: pulling functions from sourced files
# v2.0 - updated for V4.13 of pipeline
# v1.3 - Options for plotting by sample
# v1.1 - Outputs QC stats for singlets; single cell heatmap uses cap + 2

###
myPath <- .libPaths()
myPath <- c(myPath,'/mnt/bin/R/library/')
.libPaths(myPath)
source('/home/snagaraja/comp_tools/scRNA_scATAC_scripts/function_code/singlet_clustering.shared.v1.6.R')
source('/home/snagaraja/comp_tools/scRNA_scATAC_scripts/function_code/singlet_clustering.RNA.v1.4.R')

#############################################################
###################### 01 - DECLUMPING ######################
#############################################################

umi <- readRDS(umiFile)
cat('Mean transcripts per cell:',mean(Matrix::colSums(umi)),'\n')
umi_filt <- umi[,(Matrix::colSums(umi) > umi_thresh)]
cat('Number of cells passing UMI threshold',paste0('(',umi_thresh,')'),':',ncol(umi_filt),'out of',ncol(umi),'\n')

qc <- read.csv(qc_file,header=T)
colnames(qc)[grepl('_libsize',colnames(qc))] <- 'libsize'
qc$barcode <- with(qc,paste(R1,R2,R3,P5,sep=','))
if (sum(!(colnames(umi_filt) %in% qc$barcode)) > 0){
  cat('There are cells in the UMI matrix not found in QC files!\n')
}
qc_sub1 <- qc[(qc$barcode %in% colnames(umi_filt)),]
qc_sub1 <- qc_sub1[match(colnames(umi_filt),qc_sub1$barcode),]

if (two_thresh){
  keep.cells <- splitThreshold( counts = umi_filt,
                                sizes = qc_sub1$libsize,
                                thresh1 = umi_thresh1,
                                thresh2 = umi_thresh2,
                                marker = thresh_marker
                                )
  dev.off()
  cat(sum(keep.cells),'cells of',length(keep.cells),'total retained after split thresholding\n')
  umi_filt <- umi_filt[,keep.cells]
  qc_sub1 <- qc_sub1[keep.cells,]
}

norm1 <- t(t(umi_filt)/(Matrix::colSums(umi_filt))) # normalized to total transcripts per cell
exp <- as(norm1*mean(Matrix::colSums(umi_filt)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(exp,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.all.rds'))

# filter for variant / expressed genes
exp.sd <- rowSds(exp)
names(exp.sd) <- rownames(exp)
var <- exp.sd[order(exp.sd,decreasing=T)]
l2e <- log2(t(exp[names(var)[1:5000],])+1)

# PCA
PCA <- cachePCA(	cachePath = './pcaCache',
					dataSet = as.matrix(l2e),
					center = F, scale = F)
pc.scores <- data.frame(PCA$x[,1:50])
saveRDS(pc.scores,paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.all.rds'))
saveRDS(data.frame(PCA$rotation[,1:50]),paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_rotations.all.rds'))

# top 20 PC UMAP
UMAP <- cacheUMAP(	cachePath = './umapCache',
					dataSet = pc.scores[,1:20],
					seed = 123)
umap.p <- data.frame(UMAP$layout)
colnames(umap.p) <- c('UMAP1','UMAP2')
saveRDS(umap.p,paste0(output,'.min_umi_',umi_thresh,'.pc_20.umap_coord.all.rds'))
rm(umi,PCA,UMAP)
gc(verbose=FALSE)

## Finding KNN
knn_data <- get.knn(pc.scores[,1:20],k=20)
knn <- knn_data$nn.index
write.table(knn,paste0(output,'.pc20_knn_20.all.csv'),col.names=F,row.names=F,quote=F,sep=',')

# QC filtering
qc_sub2 <- qc_sub1[(qc_sub1$barcode %in% rownames(pc.scores)),]
qc_sub2 <- qc_sub2[match(rownames(pc.scores),qc_sub2$barcode),]
qc_sub2$log10UMI <- log10(Matrix::colSums(umi_filt))

# identifying clumps
pdf(paste0(output,'.umap_all.log10UMI.pdf'),width=4.5,height=4)
g <- ggplot(cbind(qc_sub2,umap.p),aes(x=UMAP1,y=UMAP2,color=log10UMI)) + geom_point(size=umap_pt_size,stroke=0) + 
	scale_color_gradient(low='navy',high='yellow') +
		theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()
clump.ind <- getShareClumps(qc_sub2$libsize,knn,umap.p,ptSize=1)
writeLines(rownames(umap.p)[clump.ind],paste0(output,'.RNA_clump_barcodes.txt'))
writeLines(rownames(umap.p)[!(clump.ind)],paste0(output,'.RNA_declump_barcodes.txt'))
umi_declump <- umi_filt[,!(clump.ind)]
saveRDS(umi_declump,paste0(output,'.min_umi_',umi_thresh,'.raw_UMI.declump.rds'))
norm_declump <- t(t(umi_declump)/(Matrix::colSums(umi_declump))) # normalized to total transcripts per cell
exp_declump <- as(norm_declump*mean(Matrix::colSums(umi_declump)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(exp_declump,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.declump.rds'))
rm(umi_filt,exp,norm1,l2e)
gc(verbose=FALSE)

# variant genes, PCA, UMAP on declumped data
exp_declump.sd <- rowSds(exp_declump)
names(exp_declump.sd) <- rownames(exp_declump)
var_declump <- exp_declump.sd[order(exp_declump.sd,decreasing=T)]
l2e_declump <- log2(t(exp_declump[names(var_declump)[1:5000],])+1)
PCA.d <- cachePCA(	cachePath = './pcaCache',
                   dataSet = as.matrix(l2e_declump),
                   center = F, scale = F)
pc.scores.d <- data.frame(PCA.d$x[,1:50])
saveRDS(pc.scores.d,paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.declump.rds'))
saveRDS(data.frame(PCA.d$rotation[,1:50]),paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_rotations.declump.rds'))

UMAP.d <- cacheUMAP(	cachePath = './umapCache',
                     dataSet = pc.scores.d[,1:20],
                     seed = 123)
umap.d <- data.frame(UMAP.d$layout)
colnames(umap.d) <- c('UMAP1','UMAP2')
saveRDS(umap.d,paste0(output,'.min_umi_',umi_thresh,'.pc_20.umap_coord.declump.rds'))
rm(PCA.d,UMAP.d)
gc(verbose=FALSE)

#############################################################
################ 02 - DOUBLETS / COLLISIONS #################
#############################################################

# scrublet to identify doublets
mtx_file <- paste0(output,'.min_umi_',umi_thresh,'.raw_UMI.declump.mtx')
writeMM(t(umi_declump),file=mtx_file)
doublet_info <- runPyScrublet(mtx_file,doub_est)
rownames(doublet_info) <- colnames(umi_declump)
write.table(doublet_info,paste0(output,'.min_umi_',umi_thresh,'.scrublet_est_',doub_est,'.doublet_info.txt'),col.names=F,row.names=T,quote=F,sep='\t')
doub <- findDoublets(doublet_info[,1],umap.d,ptSize=1)
doub.index <- doub[[1]]
doub.cutoff <- doub[[2]]
writeLines(rownames(umap.d)[doub.index],paste0(output,'.RNA_doublets_',round(doub.cutoff,2),'.barcodes.txt'))

#############################################################
################### 03 - SINGLET ANALYSIS ###################
#############################################################

# Re-run PCA, UMAP on singlets
umi.s <- umi_declump[,!(doub.index)]
saveRDS(umi.s,paste0(output,'.min_umi_',umi_thresh,'.raw_UMI.singlet.rds'))
norm.s <- t(t(umi.s)/(Matrix::colSums(umi.s))) # normalized to total transcripts per cell
exp.s <- as(norm.s*mean(Matrix::colSums(umi.s)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(exp.s,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.singlet.rds'))
rm(umi_declump,norm_declump,exp_declump,l2e_declump,norm.s)
gc(verbose=FALSE)

exp.s.sd <- rowSds(exp.s)
names(exp.s.sd) <- rownames(exp.s)
var.s <- exp.s.sd[order(exp.s.sd,decreasing=T)]
l2e.s <- log2(t(exp.s[names(var.s)[1:5000],])+1)

PCA.s <- cachePCA(	cachePath = './pcaCache',
                 dataSet = as.matrix(l2e.s),
                 center = F, scale = F)
pc.scores.s <- data.frame(PCA.s$x[,1:50])
saveRDS(pc.scores.s,paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.singlet.rds'))
saveRDS(data.frame(PCA.s$rotation[,1:50]),paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_rotations.singlet.rds'))
rm(PCA.s,l2e.s)
gc(verbose=FALSE)

qc_sub.s <- qc_sub2[(qc_sub2$barcode %in% rownames(pc.scores.s)),]
qc_sub.s <- qc_sub.s[match(rownames(pc.scores.s),qc_sub.s$barcode),]
colnames(qc_sub.s)[grepl('_unique',colnames(qc_sub.s))] <- 'unique'
## plot singlet QC stats
q1 <- ggplot(qc_sub.s,aes(x=1,y=log10(libsize)))+geom_violin(fill='dodgerblue')+xlab('')
detected_genes <- Matrix::colSums(umi.s)
detected_genes <- detected_genes[match(qc_sub.s$barcode,names(detected_genes))]
qc_sub.s$det_genes <- detected_genes
q2 <- ggplot(qc_sub.s,aes(x=1,y=log10(det_genes)))+geom_violin(fill='red')+xlab('')
q3 <- ggplot(qc_sub.s,aes(x=unique,y=det_genes))+geom_point() + xlab('UMIs detected')+ylab('Genes detected')
clump_stats <- data.frame(n_cells=c(sum(clump.ind),sum(doub.index),ncol(umi.s)),type=c('clumps','doublets','singlets'))
clump_stats$frac_cells <- clump_stats$n_cells/sum(clump_stats$n_cells)
q4 <- ggplot(clump_stats,aes(y=n_cells,x=1,fill=type))+geom_bar(stat='identity',position='stack')+coord_flip()+xlim(0,2)+xlab('')+ylab('Number of cells')
q5 <- ggplot(clump_stats,aes(y=frac_cells,x=1,fill=type))+geom_bar(stat='identity',position='stack')+coord_flip()+xlim(0,2)+xlab('')+ylab('Fraction of cells')
pdf(paste0(output,'.min_umi_',umi_thresh,'.qc_stats.singlet.pdf'),width=8,height=5)
plot_grid(q1,q2,q3,nrow=1)
with(qc_sub.s,plot(log10(qc_sub.s$libsize)[order(qc_sub.s$libsize,decreasing=T)],col='firebrick',ylab='log10(library size)',xlab='Rank',main='Singlets - RNA'))
plot_grid(q4,q5,ncol=1)
dev.off()

## if multiple samples per run, enter here as "sample" column
## e.g. 
## r1.s <- as.numeric(stri_split_fixed(qc_sub.s$R1,'.',simplify=T)[,2])
## qc_sub.s$sample <- 'blank'
## qc_sub.s$sample[(r1.s<=16)] <- 'd35.1'
## qc_sub.s$sample[(r1.s>16) & (r1.s<=48)] <- 'Control-P1'
r1.s <- as.numeric(stri_split_fixed(qc_sub.s$R1,'.',simplify=T)[,2])
qc_sub.s$sample <- 'blank'
qc_sub.s$replicate <- 'blank'
qc_sub.s$batch <- 'blank'
qc_sub.s$phase <- 'blank'

write.table(qc_sub.s,paste0(output,'.min_umi_',umi_thresh,'.qc_stats.singlet.txt'),quote=F,sep='\t',row.names=F)

### check for depth in PCs
depth.check <- pcaDepthCor(pc.scores.s,qc_sub.s$libsize)

cors <- depth.check[[1]]
usePC1 <- depth.check[[2]]

pdf(paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.depth_cor.singlet.pdf'),width=5,height=4)
plot(cors[1:20],xlab='Principal Component',ylab='Correlation with log10 library size')
lines(cors[1:20])
text(2,0.9*max(cors[1:20]),labels=paste0('Use PC1: ',usePC1))
dev.off()

if (usePC1 == 'n' | usePC1 == 'N'){
  cat('Excluding PC1 from clustering\n')
  p <- pc.scores.s[,2:21]
} else {
  p <- pc.scores.s[,1:20]
}

UMAP.s <- cacheUMAP(	cachePath = './umapCache',
                   dataSet = p,
                   seed = 123)
umap.s <- data.frame(UMAP.s$layout)
colnames(umap.s) <- c('UMAP1','UMAP2')
saveRDS(umap.s,paste0(output,'.min_umi_',umi_thresh,'.pc_20.umap_coord.singlet.rds'))
rm(UMAP.s)
gc(verbose=FALSE)

# Smoothening KNN
knn_data.s <- get.knn(p,k=20)
knn.s <- knn_data.s$nn.index
write.table(knn.s,paste0(output,'.pc20_knn_20.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')

# plot depth
qc_sub.s <- cbind(qc_sub.s,umap.s)
pdf(paste0(output,'.umap_singlet.libsize.pdf'),width=4.5,height=4)
g <- ggplot((qc_sub.s),aes(x=UMAP1,y=UMAP2,color=log10(libsize))) + geom_point(size=umap_pt_size,stroke=0) + 
  scale_color_gradient(low='navy',high='yellow') +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()

# plot markers
exp.s.smooth <- smoothCellsNN(E = exp.s, K = knn.s, cores = cores)
saveRDS(exp.s.smooth,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.knn_20_smooth.singlet.rds'))
gc(verbose=FALSE)

if (!dir.exists('./umap_gene_plots')){
  dir.create('./umap_gene_plots')
}
setwd('./umap_gene_plots')
for (marker in markerList){
  pdf(paste0(output,'.umap_singlet.',marker,'.pdf'),width=4.5,height=4)
  plotUMIgene(umap.s,exp.s.smooth,marker,zCap,umap_pt_size)
  dev.off()
}
setwd('../')

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
  for (S in unique(qc_sub.s$phase)){
    cat(S,'\n')
    plotGroupDensity( G = S,
                      U = umap.s,
                      groups = qc_sub.s$phase,
                      NN = knn.s,
                      ncores = cores
    )
  }
  dev.off()
}

if (sum(colnames(qc_sub.s) == 'replicate') > 0) {
  pdf(paste0(output,'.umap_singlet.by_replicate.pdf'),width=5,height=4)
  g1 <- ggplot(shuf(qc_sub.s),aes(x=UMAP1,y=UMAP2,color=replicate)) + geom_point(size=umap_pt_size,stroke=0) +
    scale_color_viridis(discrete = T, option = "D") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g1)
  dev.off()
}


# louvain clustering
ig <- graph.empty(nrow(knn.s))
edge_list <- pbmclapply(X=1:nrow(knn.s),FUN=function(x){
  return(unlist(mapply(c,rep(x,ncol(knn.s)),knn.s[x,],SIMPLIFY=F)))
},mc.cores=cores)
ig <- add_edges(ig,unlist(edge_list))
comm <- cluster_louvain(as.undirected(ig))
clust <- data.frame(sample=rownames(p),community=as.vector(membership(comm)),stringsAsFactors = F)
write.table(clust,paste(output,'pc20_singlet.louvain_clusters.txt',sep='.'),quote=F,sep='\t',row.names=F)

u <- cbind(umap.s,as.character(clust$community))
colnames(u)[3] <- 'community'

pdf(paste0(output,'.umap_singlet.louvain_clusters.pdf'),width=5,height=4)
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

## plot by sample
if (isMultiSample){
  plotBySampleRNA(qc_sub.s,clust)
}

# by subgroup
dir.create('./by_passage/')
setwd('./by_passage/')
passages <- c('P1','P4')
for (passage in passages){
  cat(passage,'\n')
  passage_ind <- grepl(passage,qc_sub.s$sample)
  qc_p <- qc_sub.s[passage_ind,]
  pc_scores_p <- p[passage_ind,]
  
  UMAP.pass <- cacheUMAP(	cachePath = './umapCache',
                       dataSet = pc_scores_p,
                       seed = 123)
  umap.pass <- data.frame(UMAP.pass$layout)
  colnames(umap.pass) <- c('UMAP1','UMAP2')
  saveRDS(umap.pass,paste0(output,'.passage_',passage,'.umap_coord.singlet.rds'))
  rm(UMAP.pass)
  gc(verbose=FALSE)
  
  # Smoothening KNN
  knn_data.pass <- get.knn(pc_scores_p,k=20)
  knn.pass <- knn_data.pass$nn.index
  write.table(knn.pass,paste0(output,'.passage_',passage,'.pc20_knn_20.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')
  
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
  
  if (sum(colnames(qc_sub.s) == 'phase') > 0){
    pdf(paste0(output,'.umap_singlet.passage_',passage,'.phase_knn_density.pdf'),width=5,height=4)
    for (S in unique(qc_p$phase)){
      cat(S,'\n')
      plotGroupDensity( G = S,
                        U = umap.pass,
                        groups = qc_p$phase,
                        NN = knn.pass,
                        ncores = cores
      )
    }
    dev.off()
  }
}

# top genes per cluster
exp.z <- t(scale(t(exp.s.smooth))) #z-score gene expression
uniq_clust <- unique(clust$community)

pdf(paste0(output,'.pc20_singlet.louvain_top_genes.pdf'),height=20,width=6)
exp.z.mean <- plotByClust(exp.z,clust,uniq_clust,nPlot = nTop)
dev.off()

write.table(exp.z.mean,paste0(output,'.pc20_singlet.louvain_expr_z.csv'),sep=',',quote=F)
top.mean <- matrix(ncol=2,nrow=0,data=0)
colnames(top.mean) <- c('cluster','gene')
for (i in 1:ncol(exp.z.mean)){
  sort <- exp.z.mean[order(exp.z.mean[,i],decreasing=T),]
  top.mean <- rbind(top.mean,data.frame(cluster=colnames(exp.z.mean)[i],gene=rownames(sort[1:nTop,])))
}
write.table(top.mean,paste0(output,'.pc20_singlet.louvain_top_genes.csv'),sep=',',row.names=F,quote=F)
cat('Done!\n')
