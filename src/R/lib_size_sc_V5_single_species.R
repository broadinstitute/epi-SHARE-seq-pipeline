#!/usr/bin/Rscript

args <- commandArgs()
elbow_knee_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)
  # Distance from point to line
  distances <- numeric(length(x_values))
  
  for(i in 1:length(x_values)) {
    distances[i] <- abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2)
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}

unfilt <- args[6]
filt <- args[7]
Cutoff <- as.integer(args[8])
Genome1 <- args[9]
libtype <- args[10]

hg_unfil <- read.csv(unfilt, header=T)
hg_fil <- read.csv(filt, header=T)

# ensure rows are in same order
hg_unfil <- hg_unfil[order(hg_unfil$R1, hg_unfil$R2, hg_unfil$R3),]
hg_fil <- hg_fil[order(hg_fil$R1, hg_fil$R2, hg_fil$R3),]

Df <- cbind(hg_fil, hg_unfil$fragments)

Genome1_unique <- paste(Genome1, "unique", sep ="_")
Genome1_total <- paste(Genome1, "total", sep ="_")
Genome1_dup <- paste(Genome1, "dup", sep ="_")

ncol = dim(Df)[2]
colnames(Df)[ncol] = Genome1_total
colnames(Df)[ncol-1] = Genome1_unique

Df$dup <- 1 - Df[ ,ncol-1]/Df[ ,ncol]
colnames(Df)[ncol+1] <- Genome1_dup

Genome1_libsize <- paste(Genome1, "libsize", sep ="_")
Df$temp <- 0
colnames(Df)[ncol+2] <- Genome1_libsize

# remove barcode combination with less than N reads
Df <- Df[Df[ ,ncol] >= Cutoff, ]

for (i in 1:nrow(Df)){
  if (Df[i,ncol] != Df[i,ncol-1]) {
    Df[i,ncol+2] <- round(uniroot(function(x) (x*(1-exp(-Df[i,ncol]/x))-Df[i,ncol-1]),
                                  lower = 0, upper = 10000000, tol = 1e-7)$root)}
  else{
    Df[i,ncol+2] <- Df[i,ncol]
  }
}
# for pair-end deviced library size by 2
if (libtype != "RNA"){
  if (libtype != "ATAC"){
    Df[ ,ncol+2] <- Df[ ,ncol+2]/2
  }
}
print(paste("Average ", Genome1, " lib size of sample is:", sep=""))
mean(Df[Df[, ncol+2] > Cutoff,ncol+2])

alllibeize <- sum(Df[ ,ncol+2])
duprate <- mean(Df[Df[, ncol+2] > 500,ncol+1])
allreads <- sum(Df[ ,ncol])
dupreads <- sum(Df[ ,ncol]) - sum(Df[ ,ncol-1])
out = data.frame(LIBRARY="Unknown",
                 UNPAIRED_READS_EXAMINED=0,
                 READ_PAIRS_EXAMINED=allreads ,
                 SECONDARY_OR_SUPPLEMENTARY_RDS=0,
                 UNMAPPED_READS=0,
                 UNPAIRED_READ_DUPLICATES=0,
                 READ_PAIR_DUPLICATES= dupreads,
                 READ_PAIR_OPTICAL_DUPLICATES=0,
                 PERCENT_DUPLICATION=duprate,
                 ESTIMATED_LIBRARY_SIZE=alllibeize)
fileout <- paste(basename(filt), 'dups.log.txt', sep=".")
write.table(out, fileout, quote=F, row.names=F, col.names=T, sep="\t")

# Plot lib size
options(scipen=999)
Size2 <-Df[, ncol+2];  Size2 <- Size2[Size2 >0]; # plot(Size)
file6 <- paste(basename(filt), 'libSize.png', sep=".")
png(file6, width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,1))
Size2 <- sort(Size2, decreasing = T)
data <- data.frame(1:length(Size2), Size2)
colnames(data) <- c("x","y")
elbow <- elbow_knee_finder(data$x, log10(data$y))
top_ranked <- factor(ifelse(data$x <= elbow[1], 1, 0))
plot(Size2,
     log="y", 
     xlab=paste0(" Barcode rank (", length(Size2)-elbow[1], " low quality cells)"), 
     ylab = "Reads per barcode", main = paste0(libtype, " Reads per barcode"), 
     col=c("dimgrey","darkblue")[top_ranked], 
     pch=16,
     ylim=c(1,100000))
abline(v=elbow[1], h=10^(elbow[2]))
text(elbow[1], 10^(elbow[2]),
     paste("(", elbow[1], ", ", 10^(elbow[2]), ")", sep=""),
     adj=c(-0.1,-0.5))

# plot top ranked barcodes
smooth_spline <- smooth.spline(x=data$x[1:elbow[1]], y=log10(data$y[1:elbow[1]]), spar=1)
second_deriv <- predict(smooth_spline, x=data$x[1:elbow[1]], deriv=2)

# if second derivative is all pos/all neg, plot is not "bumpy" so can use basic elbow/knee finder.
# second deriv vals can be noisy at beginning and end of graph; exclude first 10% and last 10% 
# of values when establishing uniformity of second deriv sign.
# else, look for zeros of second deriv plot surrounding abs min of second deriv;
# use these zeros as endpoints for elbow/knee finder function
ten_percent <- round(length(second_deriv$x)*0.1)
mid_second_deriv <- second_deriv$y[(ten_percent+1):(length(second_deriv$y)-ten_percent)] 
if (all(mid_second_deriv >= 0) | all(mid_second_deriv <= 0)){
  print("Finding elbow of top ranked barcodes")
  point <- elbow_knee_finder(x_values=data$x[1:elbow[1]],
                             y_values=log10(data$y[1:elbow[1]]))
} else {
  # find abs min
  abs_min_idx <- second_deriv$x[which.min(second_deriv$y)]
  # find last non-negative value before abs min 
  left_vect <- second_deriv$y[0:abs_min_idx]
  endpt_1_idx <- tail(which(left_vect >= 0), n=1)
  # find first non-positive value after abs min
  right_vect <- second_deriv$y[abs_min_idx:length(second_deriv$y)]
  endpt_2_idx <- abs_min_idx + which(right_vect >= 0)[1]
  
  # error case: revert to elbow finder.
  # used when second deriv curve looks all pos/all neg, but has some 
  # opposite sign values at start/end that escape mid_second_deriv cutoff
  if (length(endpt_1_idx)==0 | length(endpt_2_idx)==0){
    print("Finding elbow of top ranked barcodes")
    point <- elbow_knee_finder(x_values=data$x[1:elbow[1]],
                               y_values=log10(data$y[1:elbow[1]]))
  } else {
    # find knee between these points
    print("Finding knee of top ranked barcodes")
    point <- elbow_knee_finder(x_values=(endpt_1_idx:endpt_2_idx),
                               y_values=log10(Size2[endpt_1_idx:endpt_2_idx]))
  }
}
top_ranked <- factor(ifelse(data$x <= point[1], 1, 0))
plot(Size2[1:elbow[1]],
     log="y",
     xlab="Barcode rank",
     ylab = "Reads per barcode",
     main = paste0(libtype, " Reads per (top ranked) barcode"),
     col=c("dimgrey","darkblue")[top_ranked],
     pch=16,
     ylim=c(1,100000))
abline(v = point[1], h = 10^(point[2]))
text(point[1], 10^(point[2]),
     paste("(", point[1], ", ", 10^(point[2]), ")", sep=""),
     adj=c(-0.1,-0.5))
garbage <- dev.off()

## Removing for now as we've been told it's hard to interpret
#file12 <- paste(basename(filt), 'FilvsUnfil.pdf', sep=".")
#pdf(file12)
#plot(Df[ , ncol-1],Df[, ncol], xlab=paste(Genome1, " unique molecules", sep=""), ylab = "total reads", main = "Saturation", col="darkblue", pch=16)
#abline(1,1, col="blue")
#garbage <- dev.off()

temp <- paste(basename(filt),"libsize.counts.csv", sep=".")
write.csv(Df, temp)

# C/X = 1 - exp( -N/X )
#
# where X = number of distinct molecules in library
# N = number of read pairs
# C = number of distinct fragments observed in read pairs"
