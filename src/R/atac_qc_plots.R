#!/usr/bin/Rscript

### Takes ATAC barcode metadata tsv file, and outputs barcode rank plots as a png file.

library(ggplot2)

## Import helper functions
source("/usr/local/bin/barcode_rank_functions.R")

## Get arguments, read input
args <- commandArgs()

barcode_metadata_file <- args[6]
fragment_min_cutoff <- as.integer(args[7])
hist_min_fragment <- as.integer(args[8])
hist_max_fragment <- as.integer(args[9])
fragment_rank_plot_file <- args[10]
fragment_histogram_plot_file <- args[11]

barcode_metadata <- read.table(barcode_metadata_file, header=T)

## Get plot inputs

# Impose fragment cutoff, sort in decreasing order, assign rank
# 1 fragment = 2 reads
fragment <- barcode_metadata$unique
fragment_filtered <- fragment[fragment >= fragment_min_cutoff]
fragment_filtered_sort <- sort(fragment_filtered, decreasing=T)
fragment_rank <- 1:length(fragment_filtered_sort)

# Find elbow/knee of fragment barcode rank plot and top-ranked fragment barcode rank plot
fragment_points <- get_elbow_knee_points(x=fragment_rank, y=log10(fragment_filtered_sort))
# For each valid plot, make factor for coloring plot points
if (length(fragment_points) > 0) { # Elbow found in first plot
  fragment_plot1 <- TRUE
  is_top_ranked_fragment <- factor(ifelse(fragment_rank <= fragment_points[1], 1, 0))
  if (length(fragment_points) > 2) { # Elbow/knee found in second plot
    fragment_plot2 <- TRUE
    fragment_top_rank <- fragment_rank[1:fragment_points[1]]
    fragment_top_fragment <- fragment_filtered_sort[1:fragment_points[1]]
    is_top_top_ranked_fragment <- factor(ifelse(fragment_top_rank <= fragment_points[3], 1, 0))
  } else {
    fragment_plot2 <- FALSE
  }
} else {
    fragment_plot1 <- FALSE
}

## Generate plots

options(scipen=999)

# Make fragment barcode rank plots
png(fragment_rank_plot_file, width=8, height=8, units='in', res=300)
par(mfrow = c(2,1))

# Plot 1 (all barcodes passing fragment filter vs log10(fragments))
if (fragment_plot1) {
  plot(x=fragment_rank,
       y=fragment_filtered_sort,
       log="xy",
       xlab=paste0(" Barcode rank (", length(fragment_rank)-fragment_points[1], " low quality cells)"),
       ylab="Fragments per barcode",
       main="ATAC Fragments per Barcode",
       col=c("dimgrey","darkblue")[is_top_ranked_fragment],
       pch=16,
       ylim=c(1,500000))
  abline(v=fragment_points[1], h=10^(fragment_points[2]))
  text(fragment_points[1], 10^(fragment_points[2]),
       paste0("(", fragment_points[1], ", ", 10^(fragment_points[2]), ")"),
       adj=c(-0.1,-0.5))
}

# Plot 2 (top ranked barcodes vs log10(fragments))
if (fragment_plot2) {
  plot(x=fragment_top_rank,
       y=fragment_top_fragment,
       log="xy",
       xlab="Barcode rank",
       ylab="Fragments per barcode",
       main="ATAC Fragments per Top-Ranked Barcode",
       col=c("dimgrey","darkblue")[is_top_top_ranked_fragment],
       pch=16,
       ylim=c(1,500000))
  abline(v=fragment_points[3], h=10^(fragment_points[4]))
  text(fragment_points[3], 10^(fragment_points[4]),
       paste("(", fragment_points[3], ", ", 10^(fragment_points[4]), ")", sep=""),
       adj=c(-0.1,-0.5))
}
dev.off()

# Make fragment count histogram
png(fragment_histogram_plot_file, width=6, height=4, units='in', res=300)

ggplot(mapping=aes(fragment[fragment >= hist_min_fragment])) + 
  geom_histogram(binwidth=100) + 
  xlim(0, hist_max_fragment) +
  xlab("Fragments per barcode")

dev.off()