library(reticulate)
use_python("/usr/bin/python3")

runPyScrublet <- function(mm_file, #UMI matrix file
                          est=0.06
) {
  py_run_string("import os")
  py_run_string("import scipy.io")
  py_run_string("import numpy as np")
  py_run_string("import scrublet")
  py_run_string("import argparse")
  cat("Reading in data..\n")
  py_run_string(paste0("mat = scipy.io.mmread('",mm_file,"')"))
  cat('Running scrublet...\n')
  py_run_string(paste0("scrub = scrublet.Scrublet(mat, expected_doublet_rate=float(",est,"))"))
  py_run_string("doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)")
#  py_run_string("is_doublet = [int(d) for d in predicted_doublets]")
#  py_run_string("doublet_info = np.array([list(doublet_scores),is_doublet])")
 py_run_string("doublet_info = np.array([list(doublet_scores)])")
  return(t(py$doublet_info))
}

findDoublets <- function(scores, #vector of doublet score 
                         umap, #matched UMAP coordinate matrix, columns UMAP1 and UMAP2
                         cutoff=0.49, #default cut off
                         ptSize=0.5
){
  u <- cbind(umap,scores)
  s.max <- max(scores)
  
  g <- ggplot(u,aes(x=UMAP1,y=UMAP2,color=scores))+ geom_point(size=ptSize,stroke=0) + 
    scale_color_gradient(low='beige',high='navy') +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
    labs(color='Doublet score')
  print(g)
  
  adv <- readline(prompt='Hit enter to advance')
  dev.off()
  
  
  cat('Starting doublet score cutoff:',cutoff,'\n')
  hist(scores,xlab='Doublet score',main='RNA-seq Doublet Scores - Scrublet',breaks=40)
  abline(v=cutoff,lty=2,col='red')
  cat('Number of doublets called:',sum(scores > cutoff),'out of',length(scores),'\n')
  cont <- readline(prompt="Adjust cutoff? (y/n): ")
  dev.off()
  
  while (cont == 'y' | cont == 'Y'){
    hist(scores,xlab='Doublet score',main='RNA-seq Doublet Scores - Scrublet',breaks=40)
    abline(v=cutoff,lty=2,col='red')
    cutoff <- as.numeric(readline(prompt='Enter doublet score cutoff to use: '))
    dev.off()
    g1 <- ggplot(u,aes(x=UMAP1,y=UMAP2,color=scores))+ geom_point(size=1,stroke=0) + 
      scale_color_gradient(low='beige',high='navy',limits=c(min(scores),max(scores))) + ggtitle('Declumped') +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
      labs(color='Doublet score')
    g2 <- ggplot(u[(scores <= cutoff),],aes(x=UMAP1,y=UMAP2,color=scores[(scores <= cutoff)]))+ geom_point(size=1,stroke=0) + 
      scale_color_gradient(low='beige',high='navy',limits=c(min(scores),max(scores))) + ggtitle('Doublets out') +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
      labs(color='Doublet score')
    print(plot_grid(g1,g2,nrow=1))
    cat('Number of doublets called:',sum(scores > cutoff),'out of',length(scores),'\n')
    cont <- readline(prompt="Adjust cutoff? (y/n): ")
  }
  dev.off()
  
  pdf(paste0(output,'.doublet_umap.declump.pdf'),width=13,height=6)
  print(plot_grid(g1,g2,nrow=1))
  dev.off()
  
  pdf(paste0(output,'.doublet_hist.declump.pdf'),width=5,height=3)
  hist(scores,xlab='Doublet score',main='RNA-seq Doublet Scores - Scrublet',breaks=40)
  abline(v=cutoff,lty=2,col='red')
  dev.off()
  
  return(list((scores > cutoff),cutoff))
}

pcaDepthCor <- function(pc, #matrix of principal components, columns are PCs
                        sizes #vector library sizes, match rows of pc
){
  c <- unlist(lapply(X=1:ncol(pc),FUN=function(x){
    return(cor(log10(sizes),pc[,x]))
  }))
  plot(c[1:20],xlab='Principal Component',ylab='Correlation with log10 library size')
  lines(c[1:20])
  use <- readline(prompt='Use PC1 for UMAP and clustering? (y/n): ')
  dev.off()
  return(list(c,use))
}

plotBySampleRNA <- function(qc_mat, #metadata matrix with columns UMAP1, UMAP2, sample
                         clust_mat #clustering matrix, columns 'sample' and 'community'
){
  pdf(paste0(output,'.umap_singlet.by_sample.pdf'),width=5,height=4)
  g1 <- ggplot(shuf(qc_mat),aes(x=UMAP1,y=UMAP2,color=sample)) + geom_point(size=umap_pt_size,stroke=0) +
    scale_color_viridis(discrete = T, option = "D") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g1)
  dev.off()
  print(g1)
  adv <- readline(prompt='Hit enter to advance')
  dev.off()
  
  clust.all <- unique(clust_mat$community)
  sample.all <- as.character(unique(qc_mat$sample))
  count.list <- lapply(X=clust.all, FUN=function(x){
    clust.cells <- as.character(clust_mat[(clust_mat$community == x),1])
    counts <- table(qc_mat$sample[qc_mat$barcode %in% clust.cells])
    return(counts[match(sample.all,names(counts))])
  })
  cnt <- do.call('rbind',count.list)
  rownames(cnt) <- paste0('comm_',clust.all)
  write.table(cnt,paste0(output,'.cluster_singlet.by_sample.csv'),sep=',',quote=F)
  
  freq.list <- lapply(X=1:nrow(cnt),FUN=function(x){return(cnt[x,]/sum(cnt[x,]))})
  cnt.freq <- do.call('rbind',freq.list)
  rownames(cnt.freq) <- rownames(cnt)
  pdf(paste0(output,'.cluster_singlet.by_sample.pdf'),width=8,height=4)
  barplot(t(cnt.freq),legend=T,cex.names=0.5,ylab='Fraction of cells',main='Cluster Composition - By Sample')
  dev.off()
  barplot(t(cnt.freq),legend=T,cex.names=0.5,ylab='Fraction of cells',main='Cluster Composition - By Sample')
  adv <- readline(prompt='Hit enter to advance')
  
  pdf(paste0(output,'.qc_stats.singlet.by_sample.pdf'),width=6,height=4)
  qs1 <- ggplot(qc_mat,aes(x=sample,y=log10(libsize),fill=sample)) + geom_violin(scale='area',trim = F) + 
    scale_color_viridis(discrete = T, option = "D") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  qs2 <- ggplot(qc_mat,aes(x=sample,y=log10(det_genes),fill=sample)) + geom_violin(scale='area',trim = F) + 
    scale_color_viridis(discrete = T, option = "D") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(qs1)
  print(qs2)
  barplot(table(qc_mat$sample),ylab='Number of cells',xlab='sample')
  dev.off()
  
  print(qs1)
  adv <- readline(prompt='Hit enter to advance')
  print(qs2)
  adv <- readline(prompt='Hit enter to advance')
  barplot(table(qc_mat$sample),ylab='Number of cells',xlab='sample')
  adv <- readline(prompt='Hit enter to advance')
  dev.off()
}
