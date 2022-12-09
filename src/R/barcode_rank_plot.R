#!/usr/bin/Rscript

### Takes in vector of UMIs per barcode for RNA, or fragments per barcode for ATAC --
### sorted in descending order.
### Generates a png with two barcode rank plots: one with all barcodes, 
### and one with only top-ranked barcodes (determined by elbow of first plot).

## Get arguments, read input

args <- commandArgs()

input_file <- args[6]
cutoff <- as.integer(args[7])
genome <- args[8]
assay <- args[9]
output_file <- args[10]

vector <- scan(input_file)

## Define functions needed for plotting 

# Helper function to get vectors on which to call the elbow_knee_finder. 
# Takes in xy values of the curve, outputs appropriate xy vectors to be passed to elbow_knee_finder.
#
# Function computes the second derivative of the curve, and uses the shape of the second
# derivative curve to determine whether the curve has multiple "joints" (i.e. want to find knee). 
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
# whereas advanced mode is used when it is not known whether the plot needs to find an 
# elbow or a knee. 
elbow_knee_finder <- function(x, y, mode="basic") {
  # Gith advanced mode, use helper function to determine which vectors to perform calculation on
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

## Generate plot inputs

if (assay == "RNA"){
  metric <- "UMIs"
} else {
  metric <- "Fragments"
}

# Only plot barcodes that pass cutoff
filtered_vector <- vector[vector >= cutoff]
rank <- 1:length(filtered_vector)

# Inputs for plot 1 (barcode rank plot of all barcodes passing filter)
# Find elbow of plot
elbow <- elbow_knee_finder(x=rank, y=log10(filtered_vector), mode="basic")
# Make factor to color points based on if they're before or after the elbow
top_ranked <- factor(ifelse(rank <= elbow[1], 1, 0))

# Inputs for plot 2 (barcode rank plot of top-ranked barcodes)
# Get xy coordinates of top-ranked barcodes
top_ranked_x <- rank[1:elbow[1]]
top_ranked_y <- filtered_vector[1:elbow[1]]
# Find elbow or knee of top-ranked barcodes
point <- elbow_knee_finder(x=top_ranked_x, y=log10(top_ranked_y), mode="advanced")
# Make factor to color points based on if they're before or after the elbow/knee
top_top_ranked <- factor(ifelse(top_ranked_x <= point[1], 1, 0))

## Generate plots

options(scipen=999)
png(output_file, width=8, height=8, units='in', res=300)
par(mfrow = c(2,1))

# Plot 1 (barcode rank plot of all barcodes passing filter)
plot(x=rank,
     y=filtered_vector,
     log="y",
     xlab=paste0(" Barcode rank (", length(rank)-elbow[1], " low quality cells)"), 
     ylab=paste(metric, "per barcode"),
     main=paste(assay, metric, "per Barcode"), 
     col=c("dimgrey","darkblue")[top_ranked], 
     pch=16,
     ylim=c(1,100000))
abline(v=elbow[1], h=10^(elbow[2]))
text(elbow[1], 10^(elbow[2]),
     paste0("(", elbow[1], ", ", 10^(elbow[2]), ")"),
     adj=c(-0.1,-0.5))

# plot 2 (top ranked barcodes vs log10(UMIs) )
plot(x=top_ranked_x,
     y=top_ranked_y,
     log="y",
     xlab="Barcode rank",
     ylab=paste(metric, "per barcode"),
     main=paste(assay, metric, "per Top-Ranked Barcode"),
     col=c("dimgrey","darkblue")[top_top_ranked],
     pch=16,
     ylim=c(1,100000))
abline(v=point[1], h=10^(point[2]))
text(point[1], 10^(point[2]),
     paste("(", point[1], ", ", 10^(point[2]), ")", sep=""),
     adj=c(-0.1,-0.5))

dev.off()
