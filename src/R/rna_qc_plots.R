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
umi_filtered <- barcode_metadata$total_counts[barcode_metadata$total_counts >= umi_cutoff]
umi_filtered_sort <- sort(umi_filtered, decreasing=T)
umi_rank <- 1:length(umi_filtered_sort)

# Find elbow/knee of UMI barcode rank plot and top-ranked UMI barcode rank plot
umi_points <- get_elbow_knee_points(x=umi_rank, y=log10(umi_filtered_sort))
# For each valid plot, make factor for coloring plot points
if (length(umi_points) > 0) { # Elbow found in first plot
  umi_plot1 <- TRUE
  is_top_ranked_umi <- factor(ifelse(umi_rank <= umi_points[1], 1, 0))
  if (length(umi_points) > 2) { # Elbow/knee found in second plot
    umi_plot2 <- TRUE
    umi_top_rank <- umi_rank[1:umi_points[1]]
    umi_top_umi <- umi_filtered_sort[1:umi_points[1]]
    is_top_top_ranked_umi <- factor(ifelse(umi_top_rank <= umi_points[3], 1, 0))
  } else {
    umi_plot2 <- FALSE
  }
} else {
  umi_plot1 <- FALSE
}

# Impose gene cutoff, sort in decreasing order, assign rank
gene_filtered <- barcode_metadata$genes[barcode_metadata$genes >= gene_cutoff]
gene_filtered_sort <- sort(gene_filtered, decreasing=T)
gene_rank <- 1:length(gene_filtered_sort)

# Find elbow/knee of gene barcode rank plot and top-ranked gene barcode rank plot
gene_points <- get_elbow_knee_points(x=gene_rank, y=log10(gene_filtered_sort))
# For each valid plot, make factor for coloring plot points
if (length(gene_points) > 0) { # Elbow found in first plot
  gene_plot1 <- TRUE
  is_top_ranked_gene <- factor(ifelse(gene_rank <= gene_points[1], 1, 0))
  if (length(gene_points) > 2) { # Elbow/knee found in second plot
    gene_plot2 <- TRUE
    gene_top_rank <- gene_rank[1:gene_points[1]]
    gene_top_gene <- gene_filtered_sort[1:gene_points[1]]
    is_top_top_ranked_gene <- factor(ifelse(gene_top_rank <= gene_points[3], 1, 0))
  } else {
    gene_plot2 <- FALSE
  }
} else {
  gene_plot1 <- FALSE
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
       xlab=paste0(" Barcode rank (", length(umi_rank)-umi_points[1], " low quality cells)"), 
       ylab="log10(UMIs)",
       main="RNA UMIs per Barcode", 
       col=c("dimgrey","darkblue")[is_top_ranked_umi], 
       pch=16,
       ylim=c(1,100000))
  abline(v=umi_points[1], h=10^(umi_points[2]))
  text(umi_points[1], 10^(umi_points[2]),
       paste0("(", umi_points[1], ", ", 10^(umi_points[2]), ")"),
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
  abline(v=umi_points[3], h=10^(umi_points[4]))
  text(umi_points[3], 10^(umi_points[4]),
       paste("(", umi_points[3], ", ", 10^(umi_points[4]), ")", sep=""),
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
       xlab=paste0(" Barcode rank (", length(gene_rank)-gene_points[1], " low quality cells)"), 
       ylab="log10(genes)",
       main="RNA Genes per Barcode", 
       col=c("dimgrey","darkblue")[is_top_ranked_gene], 
       pch=16,
       ylim=c(1,10000))
  abline(v=gene_points[1], h=10^(gene_points[2]))
  text(gene_points[1], 10^(gene_points[2]),
       paste0("(", gene_points[1], ", ", 10^(gene_points[2]), ")"),
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
  abline(v=gene_points[3], h=10^(gene_points[4]))
  text(gene_points[3], 10^(gene_points[4]),
       paste("(", gene_points[3], ", ", 10^(gene_points[4]), ")", sep=""),
       adj=c(-0.1,-0.5))
}
dev.off()

# Make genes vs UMIs scatter plot
png(gene_umi_plot_file, width=8, height=8, units='in', res=300)

plot(x=barcode_metadata$total_counts,
     y=barcode_metadata$genes,
     xlab="UMIs",
     ylab="Genes",
     main="RNA Genes vs UMIs",
     col="darkblue",
     pch=16)

dev.off()
