#!/usr/bin/Rscript

## Define functions needed for plotting barcode rank

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
    endpt_2_idx <- abs_min_idx + which(right_vect >= 0)[1] - 1
    
    # Error cases: revert to elbow finder
    # Used when second derivative curve has both positive and negative values, 
    # but doesn't match positive-negative-positive shape expected of a knee's second derivative
    if (length(endpt_1_idx)==0 | length(endpt_2_idx)==0){
      print("Returning original vectors")
      return(list(x,y))
    } else if (is.na(endpt_1_idx) | is.na(endpt_2_idx)){
      print("Returning original vectors")
      return(list(x,y))
    } else {
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
    # smooth.spline() function used in get_vectors() requires at least 4 unique
    # x values; preempt this error
    if (length(unique(x)) < 4) {
      return(NULL)
    } else {
      xy_vects <- get_vectors(x, y)
      x <- xy_vects[[1]]
      y <- xy_vects[[2]]
    }
  }
  # Error case: return null if vectors have length 0
  if (length(x)==0 | length(y)==0) {
    return(NULL)
  }
  # Get endpoints (point with smallest x value, point with largest x value)
  endpts_df <- data.frame(x_coords=c(x[1], x[length(x)]),
                          y_coords=c(y[1], y[length(y)]))
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

# Function to find the elbow/knee of a plot, and the elbow/knee of the points 
# before the first elbow/knee (i.e. elbow/knee of all barcodes, and elbow/knee
# of top-ranked barcodes).
# Takes in xy coordinates of the plot and returns vector of four coordinates:
# xy coordinates of first elbow/knee, and xy coordinates of second elbow/knee.
get_elbow_knee_points <- function(x, y) {
  point_1 <- elbow_knee_finder(x, y, mode="basic")
  if (!is.null(point_1)) {
    point_2 <- elbow_knee_finder(x[1:point_1[1]], y[1:point_1[1]], mode="advanced")
  }
  return(c(point_1, point_2))
}
