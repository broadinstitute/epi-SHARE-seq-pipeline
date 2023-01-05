#!/usr/bin/Rscript

### Takes RNA barcode metadata tsv file, and outputs QC plots as png files.
### QC plots include barcode rank by number of UMIs (all barcodes and top-ranked barcodes),
### barcode rank by number of genes (all barcodes and top-ranked barcodes),
### and genes vs UMIs scatter plot.

## Get arguments, read input

args <- commandArgs()

barcode_metadata_file <- args[6]
umi_cutoff <- as.integer(args[7])
gene_cutoff <- as.integer(args[8])
umi_rank_plot_file <- args[9]
gene_rank_plot_file <- args[10]
gene_umi_plot_file <- args[11]

barcode_metadata <- read.table(barcode_metadata_file, header=T)

## Define functions needed for plotting 

# Helper function to get vectors on which to call the elbow_knee_finder. 
# Takes in xy values of the curve, outputs appropriate xy vectors to be passed to elbow_knee_finder.
#
# Function computes the second derivative of the curve, and uses the shape of the second
# derivative curve to determine whether the curve has multiple "joints" (i.e. if knee should be found). 
# If the second derivative is uniformly positive or uniformly negative, the curve has a single "joint", 
# and so elbow_knee_finder can be called on the original input vectors.
# Otherwise (multiple "joints"), find the zeroes of the second derivative to the left and right of the 
# absolute minimum of the second derivative.
# These will be the endpoints of the elbow_knee_finder, so return the slices of the xy vectors
# between these zeroes. 
get_vectors <- function(x, y){
  smooth_spline <- smooth.spline(x, y, spar=1)
  second_deriv <- predict(smooth_spline, x, deriv=2)
  
  # Second derivative values can be noisy at beginning and end of graph; exclude first 10% and last 10% 
  # of values when establishing uniformity of second derivative sign
  ten_percent <- round(length(second_deriv$x)*0.1)
  mid_second_deriv <- second_deriv$y[(ten_percent+1):(length(second_deriv$y)-ten_percent)]
  
  if (all(mid_second_deriv >= 0) | all(mid_second_deriv <= 0)){
    print("Returning original vectors")
    return(list(x,y)) }
  else {
    # Find absolute minimum
    abs_min_idx <- second_deriv$x[which.min(second_deriv$y)]
    # Find last non-negative value before absolute minimum
    left_vect <- second_deriv$y[0:abs_min_idx]
    endpt_1_idx <- tail(which(left_vect >= 0), n=1)
    # Find first non-positive value after absolute minimum
    right_vect <- second_deriv$y[abs_min_idx:length(second_deriv$y)]
    endpt_2_idx <- abs_min_idx + which(right_vect >= 0)[1]
    
    # Error case: revert to elbow finder
    # Used when second derivative curve looks all positive/all negative, but has some 
    # opposite sign values that escape mid_second_deriv cutoff
    if (length(endpt_1_idx)==0 | length(endpt_2_idx)==0){
      print("Returning original vectors")
      return(list(x,y)) }
    else {
      print("Returning sliced vectors")
      return(list(x[endpt_1_idx:endpt_2_idx], y[endpt_1_idx:endpt_2_idx]))
    }
  }
}

# Function to find the elbow or knee of a plot. 
# Takes in set of xy coordinates of the plot and mode, returns point which is farthest 
# from the line formed by the endpoints.
# Basic mode (default) is used when the plot is known to have only one "joint",
# whereas advanced mode is used when it is not known whether the function needs to find an 
# elbow or a knee. 
elbow_knee_finder <- function(x, y, mode="basic") {
  # With advanced mode, use helper function to determine which vectors to perform calculation on
  if (mode == "advanced") {
    xy_vects <- get_vectors(x, y)
    x <- xy_vects[[1]]
    y <- xy_vects[[2]]
  }
  
  # Get endpoints (point with smallest x value, point with largest x value)
  endpt_1_idx <- which.min(x)
  endpt_2_idx <- which.max(x)
  endpts_df <- data.frame(x_coords=c(x[endpt_1_idx], x[endpt_2_idx]),
                          y_coords=c(y[endpt_1_idx], y[endpt_2_idx]))
  # Fit line between endpoints
  fit <- lm(endpts_df$y_coords ~ endpts_df$x_coords)
  # For each point, get distance from line 
  distances <- numeric(length(x))
  for(i in 1:length(x)) {
    distances[i] <- abs(coef(fit)[2]*x[i] - y[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2)
  }
  
  # Get point farthest from line
  x_max_dist <- x[which.max(distances)]
  y_max_dist <- y[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}

## Get plot inputs

# Impose UMI cutoff, sort in decreasing order, assign rank
umi_filtered <- barcode_metadata$umis[barcode_metadata$umis >= umi_cutoff]
umi_filtered_sort <- sort(umi_filtered, decreasing=T)
umi_rank <- 1:length(umi_filtered_sort)

# Find elbow of UMI barcode rank plot, make factor for coloring points
umi_elbow <- elbow_knee_finder(x=umi_rank, y=log10(umi_filtered_sort), mode="basic")
is_top_ranked_umi <- factor(ifelse(umi_rank <= umi_elbow[1], 1, 0))

# Find elbow/knee of top-ranked UMI barcode rank plot, make factor for coloring points
umi_top_rank <- umi_rank[1:umi_elbow[1]]
umi_top_umi <- umi_filtered_sort[1:umi_elbow[1]]
umi_point <- elbow_knee_finder(x=umi_top_rank, y=log10(umi_top_umi), mode="advanced")
is_top_top_ranked_umi <- factor(ifelse(umi_top_rank <= umi_point[1], 1, 0))

# Impose gene cutoff, sort in decreasing order, assign rank
gene_filtered <- barcode_metadata$genes[barcode_metadata$genes >= gene_cutoff]
gene_filtered_sort <- sort(gene_filtered, decreasing=T)
gene_rank <- 1:length(gene_filtered_sort)

# Find elbow of gene barcode rank plot, make factor for coloring points
gene_elbow <- elbow_knee_finder(x=gene_rank, y=log10(gene_filtered_sort), mode="basic")
is_top_ranked_gene <- factor(ifelse(gene_rank <= gene_elbow[1], 1, 0))

# Find elbow/knee of top-ranked gene barcode rank plot, make factor for coloring points
gene_top_rank <- gene_rank[1:gene_elbow[1]]
gene_top_gene <- gene_filtered_sort[1:gene_elbow[1]]
gene_point <- elbow_knee_finder(x=gene_top_rank, y=log10(gene_top_gene), mode="advanced")
is_top_top_ranked_gene <- factor(ifelse(gene_top_rank <= gene_point[1], 1, 0))


## Generate plots

options(scipen=999)

# Make UMI barcode rank plots
png(umi_rank_plot_file, width=8, height=8, units='in', res=300)
par(mfrow = c(2,1))

# Plot 1 (all barcodes passing UMI filter vs log10(UMIs))
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

# Plot 2 (top ranked barcodes vs log10(UMIs))
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

dev.off()

# Make gene barcode rank plots
png(gene_rank_plot_file, width=8, height=8, units='in', res=300)
par(mfrow = c(2,1))

# Plot 1 (all barcodes passing gene filter vs log10(genes))
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

# Plot 2 (top ranked barcodes vs log10(genes))
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
