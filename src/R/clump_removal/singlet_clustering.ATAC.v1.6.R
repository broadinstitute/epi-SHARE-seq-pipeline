library(GenomicRanges)
library(cisTopic)
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BiocParallel)
library(dplyr)

plotAtacQC <- function( SE,
                        var1='libsize',
                        thresh1 = 1000,
                        frip.thresh = 0.2
                        ){
  if (var1 == 'libsize'){
    with(colData(SE),plot(log10(libsize),FRIP,pch=16,cex=0.1,xlab='log10(library size)'))
    pf <- with(colData(SE),(libsize>thresh1 & FRIP>frip.thresh))
  } else if (var1 == 'depth'){
    with(colData(SE),plot(log10(depth),FRIP,pch=16,cex=0.1,xlab='log10(reads in peaks)'))
    pf <- with(colData(SE),(depth>thresh1 & FRIP>frip.thresh))
  } else {
    cat('Variable 1 not recognized\n')
  }
  abline(v=log10(thresh1),lty=2,col='red')
  abline(h=frip.thresh,lty=2,col='red')
  mtext(side=3,line=-1,paste0('Cells passing QC: ',sum(pf),'/',ncol(SE)),adj=1,col='red',cex=0.75)
  
  return(pf)
}

stripAd1 <- function(bc, split.char=',',join.char=','){
  split_bc <- stri_split_fixed(bc,split.char,simplify=T)
  if (length(split_bc) > 4){
    r1 <- paste0('R1.',split_bc[grep('R1',split_bc)+1])
    r2 <- paste0('R2.',split_bc[grep('R2',split_bc)+1])
    r3 <- paste0('R3.',split_bc[grep('R3',split_bc)+1])
  } else {
    r1 <- split_bc[grep('R1',split_bc)]
    r2 <- split_bc[grep('R2',split_bc)]
    r3 <- split_bc[grep('R3',split_bc)]
  }
  return(paste(r1,r2,r3,sep=join.char))
}

matchShareSamples <- function(atac_bc,rna_bc,isMulti=FALSE,adList=NULL){
  if (isMulti){
    ## match one sample at a time
    a.samples <- c()
    r.samples <- c()
    for (i in 1:length(adList)){
      ## Convert Ad1.xx to P1.xx
      a.sub <- atac_bc[grepl(adList[[i]][2],atac_bc)]
      r.sub <- rna_bc[grepl(adList[[i]][1],rna_bc)]
      a.strip <- unlist(lapply(X=a.sub,FUN=stripAd1))
      r.strip <- unlist(lapply(X=r.sub,FUN=stripAd1))
      a.keep <- (a.strip %in% r.strip)
      r.keep <- match(a.strip[a.keep],r.strip)
      a.samples <- append(a.samples,a.sub[a.keep])
      r.samples <- append(r.samples,r.sub[r.keep])
    }
  } else {
    a.strip <- unlist(lapply(X=atac_bc,FUN=stripAd1))
    r.strip <- unlist(lapply(X=rna_bc,FUN=stripAd1))
    a.keep <- (a.strip %in% r.strip)
    r.keep <- match(a.strip[a.keep],r.strip)
    a.samples <- atac_bc[a.keep]
    r.samples <- rna_bc[r.keep]
  }
  return(list(a.samples,r.samples))
}


getCisTopics <- function(SE,pdf_suffix,nTopics=80,cores=4){
  gr <- rowRanges(SE)
  counts <- assays(SE)$counts
  rownames(counts) <- paste(seqnames(gr),ranges(gr),sep=':')
  cis <- createcisTopicObject(counts, project.name='ATAC')
  rownames(cis@cell.data) <- colnames(SE)
  cis <- addCellMetadata(cis, cell.data = colData(SE))
  if (length(nTopics) > 1){
    topic_vect <- nTopics
  } else {
    N <- ifelse(nTopics > 10,nTopics,10)
    topic_vect <- unique(append(seq(10,N,by=10),N))
  }
  
  cat('Determining number of models...\n')
  cat(' Number of cores:',cores,'\n')
  cis <- runCGSModels(cis, 
                   topic=topic_vect, 
                   seed=123, 
                   nCores=cores, 
                   burnin = 120, 
                   iterations = 150, 
                   addModels=F)
  pdf(paste0(output,'.cisTopics_model_selection.',pdf_suffix,'.pdf'),width=5,height=4)
  cis <- selectModel(cis, type='maximum')
  dev.off()
  cat('Getting regions data...\n')
  cis <- getRegionsScores(cis,method='Z-score',scaled=TRUE)
  topics <- t(modelMatSelection(cis,target='cell',method='Z-score'))
  region.data <- cis@region.data
  return(list(topics,region.data))
}

plotBySampleATAC <- function(qc_mat,clust_mat){
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
  
  pdf(paste0(output,'.qc_stats.singlet.by_sample.pdf'),width=ceiling(0.6*length(clust.all)),height=4)
  qs1 <- ggplot(qc_mat,aes(x=sample,y=log10(libsize),fill=sample)) + geom_violin(scale='area',trim = F) +
    scale_color_viridis(discrete = T, option = "D") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  qs2 <- ggplot(qc_mat,aes(x=sample,y=FRIP,fill=sample)) + geom_violin(scale='area',trim = F) +
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

plotRNAonATAC <- function( rna.mat,
                           atac.umap,
                           gene,
                           isMulti = T,
                           adList = NULL,
                           ptSize = 1
){
  matched <- matchShareSamples( atac_bc = rownames(umap.s),
                                rna_bc = colnames(rna.mat),
                                isMulti = isMulti,
                                adList = adList
  )
  u.sub <- atac.umap[matched[[1]],]
  g <- rna.mat[gene,]
  u.sub$gene <- g[matched[[2]]]
  u.sub <- u.sub[order(u.sub$gene,decreasing=F),]
  g1 <- ggplot(u.sub,aes(x=UMAP1,y=UMAP2,color=gene)) + geom_point(size=ptSize,stroke=0) +
    scale_color_gradientn(colors=jdb_palette('brewer_fire')) +
    labs(color=gene) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g1)
}

tfExprCor <- function( scores,             # chromVAR motif scores
                       expr,               # normalized UMI matrix
                       method='spearman',
			multiSample=T,
                       adList=NULL,
                       cores=1,
                       nLabel=10,
                       addGenes=NULL
                       ) {
  matchBC <- matchShareSamples( atac_bc = colnames(scores),
                                rna_bc = colnames(expr),
                                isMulti=multiSample,
                                adList = adList
				)
  cat('ATAC cells matched to RNA:',length(matchBC[[1]]),'/',ncol(scores),'\n')
  s.sub <- scores[,matchBC[[1]]]
  e.sub <- expr[,matchBC[[2]]]
  cat('Motifs with matching genes in expression  matrix:',sum(rownames(s.sub) %in% rownames(e.sub)),'/',nrow(s.sub),'\n')
  s.sub <- s.sub[rownames(s.sub) %in% rownames(e.sub),]
  cat('Correlating...\n')
  cor.list <- pbmclapply(X=rownames(s.sub),FUN=function(x){cor(s.sub[x,],e.sub[x,],method=method)},mc.cores=cores)
  hist(unlist(cor.list),xlab=paste0('correlation (',method,')'),main='TF-Expression Correlation',breaks=40)

  df <- data.frame(TF=rownames(s.sub),cor=unlist(cor.list))
  df$sd <- apply(s.sub[df$TF,],1,sd)

  top.var <- df$TF[order((df$sd)^2,decreasing=T)][1:nLabel]
  top.cor <- df$TF[order(abs(df$cor),decreasing=T)][1:nLabel]
  top.tf <- unique(c(top.var,top.cor,addGenes))
  with(df,plot(cor,sd^2,pch=16,cex=0.5,
               xlab=paste0('TF motif score to gene expression correlation (',method,')'),
               ylab='TF motif score variance'
               )
       )
  with(df[(df$TF %in% top.tf),],text(cor,sd^2,labels=TF,cex=0.75,
                                     col=ifelse(cor>0,'red','blue')
                                     )
       )      
  return(df)
}

plotTopicsUMAP <- function(U1, # UMAP coordinates for first dataset, cells x UMAP1/2
                           topicScoreList, #list of topic score matrices for each UMAP, topics x cells
                           topicList, #vector topics to plot
                           zCap=3, # cap Z-scores for visualization
                           U2=NULL, #optional second UMAP to plot on, must give second matrix in topicScoreList
                           U3=NULL, #optional third UMAP to plot on, must give third matrix in topicScoreList
                           pause=F #pause will plotting each topic in list
                           ){
  for (topic in topicList){
    u1 <- cbind(U1,topicScoreList[[1]][topic,])
    colnames(u1)[3] <- 'Topic'
    u1$Topic[u1$Topic > zCap] <- zCap
    u1$Topic[u1$Topic < (-zCap)] <- (-zCap)
    
    if (!is.null(U2)){
      u2 <- cbind(U2,topicScoreList[[2]][topic,])
      colnames(u2)[3] <- 'Topic'
      u2$Topic[u2$Topic > zCap] <- zCap
      u2$Topic[u2$Topic < (-zCap)] <- (-zCap)
    }
    
    if (!is.null(U3)){
      u3 <- cbind(U3,topicScoreList[[3]][topic,])
      colnames(u3)[3] <- 'Topic'
      u3$Topic[u3$Topic > zCap] <- zCap
      u3$Topic[u3$Topic < (-zCap)] <- (-zCap)
    }
    
    g1 <- ggplot(shuf(u1),aes(x=UMAP1,y=UMAP2,color=Topic)) + geom_point(size=1,stroke=0) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
      scale_color_gradientn(colors=jdb_palette('solar_extra')) +
      labs(color=topic)
    if (!is.null(U2)){
      g2 <- ggplot(shuf(u2),aes(x=UMAP1,y=UMAP2,color=Topic)) + geom_point(size=0.2,stroke=0) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
        scale_color_gradientn(colors=jdb_palette('solar_extra')) +
        labs(color=topic)
      if (!is.null(U3)){
        g3 <- ggplot(shuf(u3),aes(x=UMAP1,y=UMAP2,color=Topic)) + geom_point(size=0.2,stroke=0) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
          scale_color_gradientn(colors=jdb_palette('solar_extra')) +
          labs(color=topic)
        print(plot_grid(g1,g2,g3,nrow=1))
      } else {
        print(plot_grid(g1,g2,nrow=1))
      }
    } else {
      print(g1)
    }
    if (pause){
      cont <- readline(prompt = 'Press enter to continue')
    }
  }
}
