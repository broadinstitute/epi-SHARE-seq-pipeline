library(BuenRTools)
library(ggplot2)
library(cowplot)
library(viridis)
library(FNN)
library(pbmcapply)
library(igraph)
library(pheatmap)
library(stringi)
library(BuenColors)
library(matrixStats)

splitThreshold <- function(counts, # genes x cells matrix, raw UMI counts
                           sizes, # library sizes corresponding to cells
                           marker, # marker to split cells on
                           thresh1, #threshold for positive population
                           thresh2, # threshold for negative population
                           K=5, #KNN for smoothening
                           ask = T # ask for cutoff
                           ){
  
  if (length(sizes) != ncol(counts)){
    cat('Length of sizes vector does not match cells in matrix\n')
  } else {
    cat('Prepping data for PCA...\n')
    norm_all <- t(t(counts)/(Matrix::colSums(counts))) # normalized to total transcripts per cell
    exp_all <- as(norm_all*mean(Matrix::colSums(counts)),'sparseMatrix') # scale to average transcripts across all cells
    exp_all.sd <- rowSds(exp_all)
    names(exp_all.sd) <- rownames(exp_all)
    var_all <- exp_all.sd[order(exp_all.sd,decreasing=T)]
    l2e_all <- log2(t(exp_all[names(var_all)[1:5000],])+1)
    
    PCA_all <- cachePCA(	cachePath = './pcaCache',
                         dataSet = as.matrix(l2e_all),
                         center = F, scale = F)
    pc_scores_all <- data.frame(PCA_all$x[,1:50])
    
    cat('Finding KNN and smoothening marker expression...')
    knn_data_all <- get.knn(pc_scores_all[,1:20],k=K)
    knn_all <- knn_data_all$nn.index
    c_sub <- counts[marker,]
    smooth_list_all <- pbmclapply(X=1:length(c_sub),FUN=function(x){
      return(mean(c(c_sub[x],c_sub[knn_all[x,]])))
    },mc.cores=cores)
    smooth_marker <- do.call('c',smooth_list_all)
    names(smooth_marker) <- colnames(c_sub)
  }
  
  if (ask){
    hist(log2(smooth_marker+1),main=paste('Smoothened',marker,'expression'),xlab='log2(UMI counts)+1')
    marker_cutoff <- as.numeric(readline(prompt="UMI counts for cutoff?: "))
    hist(log2(smooth_marker+1),main=paste('Smoothened',marker,'expression'),xlab='log2(UMI counts)+1')
    abline(v=log2(marker_cutoff+1),col='red',lty=2)
    cont <- readline(prompt="Set new threshold? (y/n): ")
    while (cont == 'y' | cont == 'Y'){
      marker_cutoff <- as.numeric(readline(prompt="UMI counts for cutoff?: "))
      hist(log2(smooth_marker+1),main=paste('Smoothened',marker,'expression'),xlab='log2(UMI counts)+1')
      abline(v=log2(marker_cutoff+1),col='red',lty=2)
      cont <- readline(prompt="Set new threshold? (y/n): ")
    }
  } else {
    marker_cutoff <- 1
  }
  
  keep.ind <- rep(F,length(sizes))
  keep.ind[(smooth_marker >= marker_cutoff) & (sizes >= thresh1)] <- T
  keep.ind[(smooth_marker < marker_cutoff) & (sizes >= thresh2)] <- T
  return(keep.ind)
}

getShareClumps <- function(sizes, # vector of library sizes
                           knn.mat, # KNN matrix with each row being NN for each cell in libsize vector
                           u.coord, # UMAP coordinates, columns UMAP1 and UMAP2 matching libsize vector
                           n.sd=2, # auto-threshold cutoff for clumps is 2 SD
                           ptSize=0.5
                           ){
  smooth.size.list <- pbmclapply(X=1:length(sizes),FUN=function(x){
    return(mean(sizes[c(x,knn.mat[x,])]))
  },mc.cores=cores)
  smoothsize <- unlist(smooth.size.list)

  g <- ggplot(u.coord,aes(x=UMAP1,y=UMAP2,color=smoothsize)) + geom_point(size=ptSize,stroke=0) +
    scale_color_gradient(low='navy',high='yellow',name='Smooth Lib Size') +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g)
  adv <- readline(prompt='Hit enter to advance')
  dev.off()
  
  cat('Default size threshold is +',n.sd,'SD\n')
  size_thresh <- mean(smoothsize)+(n.sd)*sd(smoothsize)
  cat('Smoothened lib size threshold (mean + ',n.sd,'SD):',size_thresh,'\n')
  hist(smoothsize,breaks=40,xlab='Library Size',main='KNN Smoothened Library Size\n')
  abline(v=size_thresh,lty=2,lwd=2,col='red')
  cat('Number of clump barcodes:',sum(smoothsize>size_thresh))
  
  cont <- readline(prompt="Manually adjust size threshold? (y/n): ")
  
  while (cont == 'y' | cont == 'Y'){
    size_thresh <- as.numeric(readline(prompt='Enter size threshold to use: '))
    hist(smoothsize,breaks=40,xlab='Library Size',main='KNN Smoothened Library Size')
    abline(v=size_thresh,lty=2,lwd=2,col='red')
    cat('Number of clump barcodes:',sum(smoothsize>size_thresh))
    cont <- readline('Change size threshold? (y/n): ')
  } 
  dev.off()
  
  pdf(paste0(output,'.smooth_lib_size.declump.pdf'),width=6,height=6)
  hist(smoothsize,breaks=40,xlab='Library Size',main='KNN Smoothened Library Size')
  abline(v=size_thresh,lty=2,lwd=2,col='red')
  dev.off()
  
  
  g1 <- ggplot(u.coord,aes(x=UMAP1,y=UMAP2,color=smoothsize)) + geom_point(size=ptSize,stroke=0) + 
    scale_color_gradient(low='navy',high='yellow',name='Smooth Lib Size') + ggtitle('All Cells') + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  g2 <- ggplot(u.coord[(smoothsize <= size_thresh),],aes(x=UMAP1,y=UMAP2,color=smoothsize[(smoothsize <= size_thresh)])) + geom_point(size=ptSize,stroke=0) + 
    scale_color_gradient(low='navy',high='yellow',name='Smooth Lib Size') + ggtitle('Declumped') + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  
  print(plot_grid(g1,g2,nrow=1))
  
  adv <- readline(prompt='Hit enter to advance')
  dev.off()
  
  pdf(paste0(output,'.smooth_lib_size.umap_all.pdf'),width=13,height=6)
  print(plot_grid(g1,g2,nrow=1))
  dev.off()
  
  return(smoothsize > size_thresh)
}

smoothCellsNN <- function( E, #expression matrix, genes x cells
                           K, # KNN matrix, each row = index of KNN
                           cellThresh = 50000, #threshold for splitting up large matrices
                           nChunks = 4, #for splitting up large matrices
                           cores = 4
){
  nCells <- ncol(E)
  if (nCells > cellThresh){
    chunkSize <- ceiling(nCells/nChunks)
    cat('Large matrix detected (greater than 50k cells), smoothening in chunks of',chunkSize,'cells\n')
    starts <- seq(1,nCells,by=chunkSize)
    ends <- starts + chunkSize - 1 
    ends[length(ends)] <- nCells
    chunkList <- mapply(c,starts,ends,SIMPLIFY=F)
    E.smooth.list <- vector('list',length=nChunks)
    for (i in 1:length(chunkList)){
      cat('Chunk',i,'\n')
      chunk <- chunkList[[i]]
      smooth.list.sub <- pbmclapply(X=chunk[1]:chunk[2],FUN=function(x){
        return(Matrix::rowMeans(cbind(E[,x],E[,K[x,]])))
      },mc.cores=cores)
      E.smooth.list[[i]] <- do.call('cbind',smooth.list.sub)
      E.smooth.list[[i]] <- as(E.smooth.list[[i]],'sparseMatrix')
      rm(smooth.list.sub)
      gc(verbose=F)
    }
    cat('Combining chunks...\n')
    E.smooth <- do.call('cbind',E.smooth.list)
    colnames(E.smooth) <- colnames(E)
    return(E.smooth)
  } else {
    smooth.list <- pbmclapply(X=1:ncol(E),FUN=function(x){
      return(Matrix::rowMeans(cbind(E[,x],E[,K[x,]])))
    },mc.cores=cores)
    E.smooth <- do.call('cbind',smooth.list)
    E.smooth <- as(E.smooth,'sparseMatrix')
    colnames(E.smooth) <- colnames(E)
    return(E.smooth)
  }
}

plotUMIgene <- function(coord, # UMAP coordinates, columns UMAP1 and UMAP2
                        expr, #expression matrix with columns matching rows of UMAP matix
                        gene, #gene to plot, must be in rownames of expr
                        cap = 3, #default cap (all values > cap are set to the cap value)
                        ptSize = 0.5,
                        palette = 'brewer_spectra'
){
  if (gene %in% rownames(expr)){
    e <- scale(expr[gene,])
    e[e>cap] <- cap
    e[e<(-cap)] <- (-cap)
    c <- cbind(coord,e)
    colnames(c)[3] <- 'gene'
    g <- ggplot(shuf(c),aes(x=UMAP1,y=UMAP2,color=gene)) + geom_point(size=ptSize,stroke=0) +
      scale_color_gradientn(colors=jdb_palette(palette)) +
      labs(color=gene) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
    print(g)
  } else {
    cat(gene,' not found in expression matrix\n')
  }

}

plotByClust <- function( d, # values to plot (expression, gene scores, motifs etc)
                         comm, #community matrix with columns "sample", "commmunity"
                         clust.unique, # vector of unique cluster values
                         nPlot = 50 #number of top "values" per cluster to plot
) {
  z_list <- lapply(X=1:length(clust.unique),FUN=function(x){
    comm_cells <- comm$sample[comm$community == clust.unique[x]]
    z_mean <- apply(d[,colnames(d) %in% comm_cells],1,mean)
  })
  z_mean <- do.call('cbind',z_list)
  colnames(z_mean) <- paste0('comm_',clust.unique)
  
  top <- c()
  for (i in 1:ncol(z_mean)){
    sort <- z_mean[order(z_mean[,i],decreasing=T),]
    top <- unique(append(top,rownames(sort)[1:nPlot]))
  }
  
  d.cap <- d
  d.cap[(d.cap > (zCap+2))] <- (zCap+2)
  z_mean.cap <- z_mean
  z_mean.cap[(z_mean.cap > zCap)] <- zCap
  annot <- data.frame(comm=as.character(clust$community),row.names=clust$sample)
  pheatmap(d.cap[top,],scale='none',color=colorRampPalette(c('navy','yellow'))(n=99),show_colnames=F,annotation_col = annot,fontsize_row=2.5)
  pheatmap(z_mean.cap[top,],scale='none',color=colorRampPalette(c('navy','yellow'))(n=99),fontsize_row=2.5)
  
  return(z_mean)
}

plotByGroup <- function(U, #UMAP coordinates
                        groups, #corresponding groupings, eg cluster, sample identity etc
                        groupsToPlot, #which of those groups to show
                        title = "",
                        cols=c('black','red','blue'),
                        pt_size = 0.5
){
  N <- length(groupsToPlot)
  if (N >= 4){
    cat('Too many groups (>=4)\n')
  } else {
    colnames(U)[1:2] <- c('UMAP1','UMAP2')
    U$color <- 'gray75'
    for (i in 1:N){
      U$color[(groups == groupsToPlot[i])] <- cols[i]
    }
    with(U,plot(UMAP1,UMAP2,col=color,cex=pt_size,pch=16,main=title))
    with(U,legend(0.6*max(UMAP1),0.9*max(UMAP2),legend=groupsToPlot,col=cols[1:N],cex=0.75,pch=16))
    with(shuf(U[(groups %in% groupsToPlot),]),points(UMAP1,UMAP2,col=color,cex=pt_size,pch=16))
  }
}

plotGroupDensity <- function( G, #name of group to plot
                              U, #UMAP coordinates
                              groups, #vector of group assignments
                              NN, #KNN matrix
                              pt_size = 1,
                              palette = 'brewer_violet',
                              ncores = 4
                              ){
  density_list <- pbmclapply(X=1:nrow(NN),FUN=function(x){
    g_vect <- c(groups[x],groups[NN[x,]])
    return(sum(g_vect == G)/length(g_vect))
  },mc.cores=ncores)
  colnames(U) <- c('UMAP1','UMAP2')
  U$density <- unlist(density_list)
  
  g1 <- ggplot(shuf(U),aes(x=UMAP1,y=UMAP2,color=100*density)) + geom_point(size=pt_size,stroke=0) +
    scale_color_gradientn(colors=jdb_palette(palette)) + labs(color=paste0('% ',G)) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g1)
}
