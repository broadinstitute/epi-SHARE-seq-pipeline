#!/usr/bin/Rscript

args <- commandArgs()
# print(args)
elbow_finder <- function(x_values, y_values) {
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

# print(paste(dir, Namehg, File1, sep="/"))
hg_unfil <- read.csv(unfilt, header=T)
hg_fil <- read.csv(filt, header=T)

Df <- cbind(hg_fil, hg_unfil$fragments)

Genome1_unique <- paste(Genome1, "unique", sep ="_")
Genome1_total <- paste(Genome1, "total", sep ="_")
Genome1_dup <- paste(Genome1, "dup", sep ="_")


ncol = dim(Df)[2]
colnames(Df)[ncol] = Genome1_total
colnames(Df)[ncol-1] = Genome1_unique
#colnames(Df) <- c("R1", "R2", "R3", "P5", Genome1_unique, Genome1_total)
#head(Df)


#Df$uniqueNuc <- Df$hg_unique + Df$mm_unique
#Df$totalNuc <- Df$hg_total + Df$mm_total
Df$dup <- 1 - Df[ ,ncol-1]/Df[ ,ncol]
colnames(Df)[ncol+1] <- Genome1_dup

Genome1_libsize <- paste(Genome1, "libsize", sep ="_")
Df$temp <- 0
colnames(Df)[ncol+2] <- Genome1_libsize
#head(Df)

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
#head(Df)
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
fileout <- paste(basename(filt), 'dups.log', sep=".")
write.table(out, fileout, quote=F, row.names=F, col.names=T, sep="\t")

# PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
# Brian_ATAC.mm10.dups.log

# Plot lib size
print("Plot Libsize, Species mixing etc...")

Size2 <-Df[, ncol+2];  Size2 <- Size2[Size2 >0]; # plot(Size)
file6 <- paste(basename(filt), 'libSize.jpg', sep=".")
jpeg(file6)
par(mfrow = c(2,1))
Size2 <- sort(Size2, decreasing = T)
elbow <- elbow_finder(1:length(Size2), log10(Size2))
plot(Size2, log="y", xlab=paste0("Barcode rank (", length(Size2)-elbow[1], " low quality cells)"), ylab = "Reads per barcode", main = "Estimated Lib Sizes per Cell", col="darkblue", pch=16)
abline(v = elbow[1])
axis(1, at=elbow[1], labels=elbow[1])
plot(Size2[1:elbow[1]], log="y", xlab="Barcode rank", ylab = "Reads per barcode", main = "Estimated Lib Sizes per high quality cell", col="darkblue", pch=16)
elbow <- elbow_finder(1:elbow[1], log10(Size2[1:elbow[1]]))
abline(v = elbow[1])
abline(h = 10^(elbow[2]))

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

