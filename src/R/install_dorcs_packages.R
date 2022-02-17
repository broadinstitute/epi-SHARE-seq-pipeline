options(Ncpus = 4)

list.of.packages <- c("dplyr", "Seurat","patchwork","ggplot2","ggrepel","reshape2","circlize","networkD3","GGally","igraph","network","foreach","iterators","hdf5r","ggrastr","BiocManager","remotes","pbmcapply","doSNOW","Rmpfr", "glue","magrittr","pillar","RcppArmadillo","reticulate","rlang","yaml","rpart")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

BiocManager::install("Biostrings", update=T, ask=F)
BiocManager::install("rtracklayer", update=T, ask=F)
BiocManager::install("GenomicRanges", update=T, ask=F)
BiocManager::install("motifmatchr", update=T, ask=F)
BiocManager::install("ComplexHeatmap", update=T, ask=F)
BiocManager::install("chromVAR", update=T, ask=F)

remotes::install_github("caleblareau/BuenColors")

install.packages('IRkernel')
IRkernel::installspec()


