#!/usr/bin/Rscript

### Takes RNA barcode metadata tsv file, and outputs QC plots as png files.
### QC plots include barcode rank by number of UMIs (all barcodes and top-ranked barcodes),
### barcode rank by number of genes (all barcodes and top-ranked barcodes),
### and genes vs UMIs scatter plot.

## Import helper functions
source("/usr/local/bin/barcode_rank_functions.R")

## Get arguments, read input
args <- commandArgs()

barcode_metadata_file <- args[6]
umi_cutoff <- as.integer(args[7])
gene_cutoff <- as.integer(args[8])
umi_rank_plot_file <- args[9]
gene_rank_plot_file <- args[10]
gene_umi_plot_file <- args[11]

barcode_metadata <- read.table(barcode_metadata_file, header=T)

## Get plot inputs

# Impose UMI cutoff, sort in decreasing order, assign rank
umi_filtered <- barcode_metadata$umis[barcode_metadata$umis >= umi_cutoff]
umi_filtered_sort <- sort(umi_filtered, decreasing=T)
umi_rank <- 1:length(umi_filtered_sort)

# Find elbow of UMI barcode rank plot
umi_elbow <- elbow_knee_finder(x=umi_rank, y=log10(umi_filtered_sort), mode="basic")
# Set flag for whether elbow was found or not 
umi_plot1 <- ifelse(is.null(umi_elbow), FALSE, TRUE)
# If elbow was found, make factor for coloring plot points and proceed to find elbow/knee
# of top-ranked UMI barcode rank plot
if (umi_plot1) {
  is_top_ranked_umi <- factor(ifelse(umi_rank <= umi_elbow[1], 1, 0))
  
  umi_top_rank <- umi_rank[1:umi_elbow[1]]
  umi_top_umi <- umi_filtered_sort[1:umi_elbow[1]]
  umi_point <- elbow_knee_finder(x=umi_top_rank, y=log10(umi_top_umi), mode="advanced")
  
  # Set flag for whether elbow/knee was found or not
  umi_plot2 <- ifelse(is.null(umi_point), FALSE, TRUE)
  # If yes, make factor for coloring plot points
  if (umi_plot2) {
    is_top_top_ranked_umi <- factor(ifelse(umi_top_rank <= umi_point[1], 1, 0))
  }
} else {
  # If elbow not found previously, don't make top-ranked UMI barcode rank plot
  umi_plot2 = FALSE
}

# Impose gene cutoff, sort in decreasing order, assign rank
gene_filtered <- barcode_metadata$genes[barcode_metadata$genes >= gene_cutoff]
gene_filtered_sort <- sort(gene_filtered, decreasing=T)
gene_rank <- 1:length(gene_filtered_sort)

# Find elbow of gene barcode rank plot, make factor for coloring points
gene_elbow <- elbow_knee_finder(x=gene_rank, y=log10(gene_filtered_sort), mode="basic")
# Set flag for whether elbow was found or not 
gene_plot1 <- ifelse(is.null(gene_elbow), FALSE, TRUE)
# If elbow was found, make factor for coloring plot points and proceed to find elbow/knee
# of top-ranked gene barcode rank plot
if (gene_plot1) {
  is_top_ranked_gene <- factor(ifelse(gene_rank <= gene_elbow[1], 1, 0))
  
  gene_top_rank <- gene_rank[1:gene_elbow[1]]
  gene_top_gene <- gene_filtered_sort[1:gene_elbow[1]]
  gene_point <- elbow_knee_finder(x=gene_top_rank, y=log10(gene_top_gene), mode="advanced")
  
  # Set flag for whether elbow/knee was found or not
  gene_plot2 <- ifelse(is.null(gene_point), FALSE, TRUE)
  # If yes, make factor for coloring plot points
  if (gene_plot2) {
    is_top_top_ranked_gene <- factor(ifelse(gene_top_rank <= gene_point[1], 1, 0))
  }
} else {
  # If elbow not found previously, don't make top-ranked UMI barcode rank plot
  gene_plot2 = FALSE
}


## Generate plots

options(scipen=999)

# Make UMI barcode rank plots
png(umi_rank_plot_file, width=8, height=8, units='in', res=300)
par(mfrow = c(2,1))

# Plot 1 (all barcodes passing UMI filter vs log10(UMIs))
if (umi_plot1) {
  plot(x=umi_rank,
       y=umi_filtered_sort,
       log="y",
       xlab=paste0(" Barcode rank (", length(umi_rank)-umi_elbow[1], " low quality cells)"), 
       ylab="log10(UMIs)",
       main="RNA UMIs per Barcode", 
       col=c("dimgrey","darkblue")[is_top_ranked_umi], 
       pch=16,
       ylim=c(1,100000))
  abline(v=umi_elbow[1], h=10^(umi_elbow[2]))
  text(umi_elbow[1], 10^(umi_elbow[2]),
       paste0("(", umi_elbow[1], ", ", 10^(umi_elbow[2]), ")"),
       adj=c(-0.1,-0.5))
}

# Plot 2 (top ranked barcodes vs log10(UMIs))
if (umi_plot2) {
  plot(x=umi_top_rank,
       y=umi_top_umi,
       log="y",
       xlab="Barcode rank",
       ylab="log10(UMIs)",
       main="RNA UMIs per Top-Ranked Barcode",
       col=c("dimgrey","darkblue")[is_top_top_ranked_umi],
       pch=16,
       ylim=c(1,100000))
  abline(v=umi_point[1], h=10^(umi_point[2]))
  text(umi_point[1], 10^(umi_point[2]),
       paste("(", umi_point[1], ", ", 10^(umi_point[2]), ")", sep=""),
       adj=c(-0.1,-0.5))
}
dev.off()


# Make gene barcode rank plots
png(gene_rank_plot_file, width=8, height=8, units='in', res=300)
par(mfrow = c(2,1))

# Plot 1 (all barcodes passing gene filter vs log10(genes))
if (gene_plot1) {
  plot(x=gene_rank,
       y=gene_filtered_sort,
       log="y",
       xlab=paste0(" Barcode rank (", length(gene_rank)-gene_elbow[1], " low quality cells)"), 
       ylab="log10(genes)",
       main="RNA Genes per Barcode", 
       col=c("dimgrey","darkblue")[is_top_ranked_gene], 
       pch=16,
       ylim=c(1,10000))
  abline(v=gene_elbow[1], h=10^(gene_elbow[2]))
  text(gene_elbow[1], 10^(gene_elbow[2]),
       paste0("(", gene_elbow[1], ", ", 10^(gene_elbow[2]), ")"),
       adj=c(-0.1,-0.5))
}

# Plot 2 (top ranked barcodes vs log10(genes))
if (gene_plot2) {
  plot(x=gene_top_rank,
       y=gene_top_gene,
       log="y",
       xlab="Barcode rank",
       ylab="log10(genes)",
       main="RNA Genes per Top-Ranked Barcode",
       col=c("dimgrey","darkblue")[is_top_top_ranked_gene],
       pch=16,
       ylim=c(1,10000))
  abline(v=gene_point[1], h=10^(gene_point[2]))
  text(gene_point[1], 10^(gene_point[2]),
       paste("(", gene_point[1], ", ", 10^(gene_point[2]), ")", sep=""),
       adj=c(-0.1,-0.5))
}
dev.off()

# Make genes vs UMIs scatter plot
png(gene_umi_plot_file, width=8, height=8, units='in', res=300)

plot(x=barcode_metadata$umis,
     y=barcode_metadata$genes,
     xlab="UMIs",
     ylab="Genes",
     main="RNA Genes vs UMIs",
     col="darkblue",
     pch=16)

dev.off()
