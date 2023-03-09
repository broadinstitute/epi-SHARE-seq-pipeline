#!/usr/bin/Rscript

### Takes ATAC barcode metadata tsv file, and outputs barcode rank plots as a png file.

## Import helper functions
source("/usr/local/bin/barcode_rank_functions.R")

## Get arguments, read input
args <- commandArgs()

barcode_metadata_file <- args[6]
fragment_cutoff <- as.integer(args[7])
fragment_rank_plot_file <- args[8]

barcode_metadata <- read.table(barcode_metadata_file, header=T)

## Get plot inputs

# Impose fragment cutoff, sort in decreasing order, assign rank
fragment_filtered <- barcode_metadata$reads_unique[barcode_metadata$reads_unique >= fragment_cutoff]
fragment_filtered_sort <- sort(fragment_filtered, decreasing=T)
fragment_rank <- 1:length(fragment_filtered_sort)

# Find elbow of fragment barcode rank plot
fragment_elbow <- elbow_knee_finder(x=fragment_rank, y=log10(fragment_filtered_sort), mode="basic")
# Set flag for whether elbow was found or not 
fragment_plot1 <- ifelse(is.null(fragment_elbow), FALSE, TRUE)
# If elbow was found, make factor for coloring plot points and proceed to find elbow/knee
# of top-ranked fragment barcode rank plot
if (fragment_plot1) {
  is_top_ranked_fragment <- factor(ifelse(fragment_rank <= fragment_elbow[1], 1, 0))
  
  fragment_top_rank <- fragment_rank[1:fragment_elbow[1]]
  fragment_top_fragment <- fragment_filtered_sort[1:fragment_elbow[1]]
  fragment_point <- elbow_knee_finder(x=fragment_top_rank, y=log10(fragment_top_fragment), mode="advanced")
  
  # Set flag for whether elbow/knee was found or not
  fragment_plot2 <- ifelse(is.null(fragment_point), FALSE, TRUE)
  # If yes, make factor for coloring plot points
  if (fragment_plot2) {
    is_top_top_ranked_fragment <- factor(ifelse(fragment_top_rank <= fragment_point[1], 1, 0))
  }
} else {
  # If elbow not found previously, don't make top-ranked fragment barcode rank plot
  fragment_plot2 = FALSE
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
       log="y",
       xlab=paste0(" Barcode rank (", length(fragment_rank)-fragment_elbow[1], " low quality cells)"), 
       ylab="log10(fragments)",
       main="ATAC Fragments per Barcode", 
       col=c("dimgrey","darkblue")[is_top_ranked_fragment], 
       pch=16,
       ylim=c(1,100000))
  abline(v=fragment_elbow[1], h=10^(fragment_elbow[2]))
  text(fragment_elbow[1], 10^(fragment_elbow[2]),
       paste0("(", fragment_elbow[1], ", ", 10^(fragment_elbow[2]), ")"),
       adj=c(-0.1,-0.5))
}

# Plot 2 (top ranked barcodes vs log10(fragments))
if (fragment_plot2) {
  plot(x=fragment_top_rank,
       y=fragment_top_fragment,
       log="y",
       xlab="Barcode rank",
       ylab="log10(fragments)",
       main="ATAC Fragments per Top-Ranked Barcode",
       col=c("dimgrey","darkblue")[is_top_top_ranked_fragment],
       pch=16,
       ylim=c(1,100000))
  abline(v=fragment_point[1], h=10^(fragment_point[2]))
  text(fragment_point[1], 10^(fragment_point[2]),
       paste("(", fragment_point[1], ", ", 10^(fragment_point[2]), ")", sep=""),
       adj=c(-0.1,-0.5))
}
dev.off()

