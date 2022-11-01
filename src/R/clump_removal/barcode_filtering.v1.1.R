library(data.table)
library(Matrix)
library(stringi)
library(pheatmap)

ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    message('CellRanger version 3+ format H5')
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    message('CellRanger version 2 format H5')
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

wellReadPlot <- function(bcList,frags,Nround){
  bcSub <- gsub(',','.',bcList)
	bcNum <- data.frame(stri_split_fixed(bcSub,'.',simplify=T)[,seq(2,2*Nround,by=2)],stringsAsFactors=F)
	for (i in 1:ncol(bcNum)){
		bcNum[,i] <- as.integer(bcNum[,i])
	}
	f <- frags[(frags$bcID %in% bcList),]
	f <- f[match(bcList,f$bcID),]
	
	if (max(bcNum) > 96){
	  rd1 <- matrix(nrow=8,ncol=24,data=0)
	  rd2 <- matrix(nrow=8,ncol=24,data=0)
	  rd3 <- matrix(nrow=8,ncol=24,data=0)
	} else {
	  rd1 <- matrix(nrow=8,ncol=12,data=0)
	  rd2 <- matrix(nrow=8,ncol=12,data=0)
	  rd3 <- matrix(nrow=8,ncol=12,data=0) 
	}
	
	for (i in 1:nrow(bcNum)){
		col1 <- ceiling(bcNum[i,1]/8)
		col2 <- ceiling(bcNum[i,2]/8)
		col3 <- ceiling(bcNum[i,3]/8)
		row1 <- ifelse(bcNum[i,1]%%8==0,8,(bcNum[i,1]%%8))
		row2 <- ifelse(bcNum[i,2]%%8==0,8,(bcNum[i,2]%%8))
		row3 <- ifelse(bcNum[i,3]%%8==0,8,(bcNum[i,3]%%8))
		
		c <- f$fragments[(f$bcID == bcList[i])]
		
		rd1[row1,col1] <- rd1[row1,col1] + c
		rd2[row2,col2] <- rd2[row2,col2] + c
		rd3[row3,col3] <- rd3[row3,col3] + c
	}
	return(list(rd1,rd2,rd3))
}

plotShareWells <- function(mat,rd){
  rowlab <- c('A','B','C','D','E','F','G','H')
  if(ncol(mat) > 12){ gap.vect <-  seq(12,12*(floor(ncol(mat)/12)-1),by=12) } else {gap.vect <- NULL}
  if (rd == 1){
    col.palette <- colorRampPalette(c('white','navy'))(n=99)
  } else if (rd == 2){
    col.palette <- colorRampPalette(c('white','firebrick'))(n=99)
  } else if (rd == 3){
    col.palette <- colorRampPalette(c('white','black'))(n=99)
  }
  pheatmap(mat,
           cluster_cols= F,
           cluster_rows= F,
           color=col.palette,
           labels_col=c(1:ncol(mat)),
           gaps_col = gap.vect,
           labels_row=rowlab,
           main=paste('Round',rd)
  )
}

getGtfGenes <- function(gtfDir,g,refgene='gencode'){
  if (sum(grepl(paste0(g,'.',refgene,'.gtf'),list.files(gtfDir))) > 0){
    cat('Reading in GTF file...\n')
    gtf <- fread(paste0(gtfDir,'/',g,'.',refgene,'.gtf'),sep='\t')
    info <- stri_split_fixed((gsub('\"','',gtf$V9)),'; ',simplify=T)
    gene_list <- pbmclapply(X=1:nrow(info),FUN=function(x){
      v <- c(info[x,])
      if (sum(grepl('gene_name',v)) == 1){
        return(stri_split_fixed(v[grepl('gene_name',v)],' ',simplify=T)[2])
      } else {
        return(NULL)
      }
    })
    genes <- unlist(unique(gene_list))
    return(list(TRUE,genes[!unlist(lapply(genes,is.null))]))
  } else {
    cat('GTF file not found for genome:',g,'\n')
    return(list(FALSE,''))
  }
}