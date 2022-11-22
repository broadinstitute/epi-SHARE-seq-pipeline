#!/usr/bin/Rscript

### Takes in filtered and unfiltered counts files, calculates library size and duplicates for
### each barcode. Outputs a counts csv with these metrics, as well as a duplicates log tsv.
### Also generates a png with two barcode rank plots: one with all barcodes, 
### and one with only top-ranked barcodes (determined by elbow of first plot).

args <- commandArgs()

unfiltered_counts_file <- args[6]
filtered_counts_file <- args[7]
cutoff <- as.integer(args[8])
genome <- args[9]
libtype <- args[10]

unfiltered_counts <- read.csv(unfiltered_counts_file, header=T)
filtered_counts <- read.csv(filtered_counts_file, header=T)

# ensure dataframes are ordered by barcode
unfiltered_counts <- unfiltered_counts[order(unfiltered_counts$R1, unfiltered_counts$R2, unfiltered_counts$R3),]
filtered_counts <- filtered_counts[order(filtered_counts$R1, filtered_counts$R2, filtered_counts$R3),]

counts_df <- cbind(filtered_counts, unfiltered_counts$fragments)

colnames(counts_df)[ncol(counts_df)-1] <- "unique"
colnames(counts_df)[ncol(counts_df)] <- "total"

## compute metrics

# remove barcode combinations that don't pass cutoff
counts_df <- counts_df[counts_df$total >= cutoff,]

# calculate duplicates
counts_df$dup <- 1 - (counts_df$unique/counts_df$total)

# calculate library size
get_libsize <- function(df){
  libsize <- rep(0, nrow(df))
  for (i in 1:nrow(df)){
    if (df$unique[i] != df$total[i]){
      libsize[i] <- round(uniroot(function(x) (x*(1-exp(-df$total[i]/x))-df$unique[i]),
                               lower = 0, upper = 10000000, tol = 1e-7)$root)}
    else{
      libsize[i] <- df$total[i]
    }
  } 
  return(libsize)
}

counts_df$libsize <- get_libsize(counts_df)

# for paired-end, divide library size by 2
if (libtype != "RNA" & libtype != "ATAC"){
  counts_df$libsize <- counts_df$libsize/2
}

print(paste("Average", genome, "lib size of sample is:", sep=" "))
mean(counts_df$libsize[counts_df$libsize > cutoff])

## generate barcode rank plots 

# function to find the elbow or knee of a plot. 
# takes in set of xy coordinates of the plot, returns point which is farthest 
# from the line formed by the endpoints.
elbow_knee_finder <- function(x, y) {
  # get endpoints (point with smallest x value, point with largest x value)
  endpt_1_idx <- which.min(x)
  endpt_2_idx <- which.max(x)
  endpts_df <- data.frame(x_coords=c(x[endpt_1_idx], x[endpt_2_idx]),
                          y_coords=c(y[endpt_1_idx], y[endpt_2_idx]))
  # fit line between endpoints
  fit <- lm(endpts_df$y_coords ~ endpts_df$x_coords)
  # for each point, get distance from line 
  distances <- numeric(length(x))
  for(i in 1:length(x)) {
    distances[i] <- abs(coef(fit)[2]*x[i] - y[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2)
  }
  
  # get point farthest from line
  x_max_dist <- x[which.max(distances)]
  y_max_dist <- y[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}

# helper function to determine shape of curve, and correspondingly what set of points
# to pass to elbow_knee_finder.
# takes in xy values and second derivative values of the curve, outputs appropriate 
# xy vectors to be passed to elbow_knee_finder.
# if the second derivative is uniformly positive or uniformly negative, the plot is not "bumpy" 
# and has a single elbow or knee, so the elbow/knee finder can be called on the original input vectors.
# otherwise ("bumpy" curve with knee), find the zeroes of the second derivative to the left and right of the 
# absolute minimum of the second derivative.
# these will be the endpoints of the elbow/knee finder, so return the slices of the xy vectors
# between these zeroes. 
elbow_knee_finder_helper <- function(x, y, second_deriv){
  # second derivative values can be noisy at beginning and end of graph; exclude first 10% and last 10% 
  # of values when establishing uniformity of second derivative sign.
  ten_percent <- round(length(second_deriv$x)*0.1)
  mid_second_deriv <- second_deriv$y[(ten_percent+1):(length(second_deriv$y)-ten_percent)]
  
  if (all(mid_second_deriv >= 0) | all(mid_second_deriv <= 0)){
    print("Returning original vectors")
    return(list(x,y)) }
  else {
    # find absolute minimum
    abs_min_idx <- second_deriv$x[which.min(second_deriv$y)]
    # find last non-negative value before absolute minimum
    left_vect <- second_deriv$y[0:abs_min_idx]
    endpt_1_idx <- tail(which(left_vect >= 0), n=1)
    # find first non-positive value after absolute minimum
    right_vect <- second_deriv$y[abs_min_idx:length(second_deriv$y)]
    endpt_2_idx <- abs_min_idx + which(right_vect >= 0)[1]
    
    # error case: revert to elbow finder.
    # used when second derivative curve looks all positive/all negative, but has some 
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

libsize <- counts_df$libsize[counts_df$libsize > 0]
libsize <- sort(libsize, decreasing=T)
libsize_df <- data.frame(1:length(libsize), libsize)
colnames(libsize_df) <- c("rank", "libsize")

# for plot 1: rank vs log10(libsize) of all barcodes
# find elbow of plot
elbow <- elbow_knee_finder(x=libsize_df$rank, y=log10(libsize_df$libsize))
# make factor to color points based on if they're before or after the elbow
top_ranked <- factor(ifelse(libsize_df$rank <= elbow[1], 1, 0))

# for plot 2: rank vs log10(libsize) of top ranked barcodes
# get xy coordinates of top ranked barcodes
top_ranked_x <- libsize_df$rank[1:elbow[1]]
top_ranked_y <- libsize_df$libsize[1:elbow[1]]
# fit spline to top ranked barcodes plot, predict second derivative
smooth_spline <- smooth.spline(x=top_ranked_x, y=log10(top_ranked_y), spar=1)
second_deriv <- predict(smooth_spline, x=top_ranked_x, deriv=2)
# find elbow/knee of the top ranked barcodes vs log10(libsize) plot
elbow_knee_vects <- elbow_knee_finder_helper(x=top_ranked_x, y=log10(top_ranked_y), second_deriv=second_deriv)
point <- elbow_knee_finder(x=elbow_knee_vects[[1]], y=elbow_knee_vects[[2]])
# make factor to color points based on if they're before or after the elbow/knee
top_top_ranked <- factor(ifelse(top_ranked_x <= point[1], 1, 0))

options(scipen=999)
libsize_plot_filename <- paste(basename(filtered_counts_file), 'libsize.png', sep=".")
png(libsize_plot_filename, width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,1))

# plot 1 (all barcodes)
plot(x=libsize_df$rank,
     y=libsize_df$libsize,
     log="y",
     xlab=paste0(" Barcode rank (", length(libsize)-elbow[1], " low quality cells)"), 
     ylab = "Reads per barcode", main = paste0(libtype, " Reads per barcode"), 
     col=c("dimgrey","darkblue")[top_ranked], 
     pch=16,
     ylim=c(1,100000))
abline(v=elbow[1], h=10^(elbow[2]))
text(elbow[1], 10^(elbow[2]),
     paste("(", elbow[1], ", ", 10^(elbow[2]), ")", sep=""),
     adj=c(-0.1,-0.5))

# plot 2 (top ranked barcodes)
plot(x=top_ranked_x,
     y=top_ranked_y,
     log="y",
     xlab="Barcode rank",
     ylab = "Reads per barcode",
     main = paste0(libtype, " Reads per (top ranked) barcode"),
     col=c("dimgrey","darkblue")[top_top_ranked],
     pch=16,
     ylim=c(1,100000))
abline(v = point[1], h = 10^(point[2]))
text(point[1], 10^(point[2]),
     paste("(", point[1], ", ", 10^(point[2]), ")", sep=""),
     adj=c(-0.1,-0.5))

dev.off()

# output duplicates log
total_libsize <- sum(counts_df$libsize)
dup_rate <- mean(counts_df$dup[counts_df$libsize > 500])
total_reads <- sum(counts_df$total)
dup_reads <- sum(counts_df$total - counts_df$unique)
dup_log_df <- data.frame(LIBRARY="Unknown",
                         UNPAIRED_READS_EXAMINED=0,
                         READ_PAIRS_EXAMINED=total_reads,
                         SECONDARY_OR_SUPPLEMENTARY_RDS=0,
                         UNMAPPED_READS=0,
                         UNPAIRED_READ_DUPLICATES=0,
                         READ_PAIR_DUPLICATES= dup_reads,
                         READ_PAIR_OPTICAL_DUPLICATES=0,
                         PERCENT_DUPLICATION=dup_rate,
                         ESTIMATED_LIBRARY_SIZE=total_libsize)
dup_log_filename <- paste(basename(filtered_counts_file), 'dups.log.txt', sep=".")
write.table(dup_log_df, dup_log_filename, quote=F, row.names=F, col.names=T, sep="\t")

# output counts df
genome_unique <- paste(genome, "unique", sep ="_")
genome_total <- paste(genome, "total", sep ="_")
genome_dup <- paste(genome, "dup", sep ="_")
genome_libsize <- paste(genome, "libsize", sep ="_")
colnames(counts_df) <- c("R1", "R2", "R3", "PKR", genome_unique, genome_total, genome_dup, genome_libsize)
counts_filename <- paste(basename(filtered_counts_file),"libsize.counts.csv", sep=".")
write.csv(counts_df, counts_filename)
