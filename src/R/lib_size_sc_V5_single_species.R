#!/usr/bin/Rscript

args <- commandArgs()
# print(args)

dir <- args[6]
Name <- args[7]
Cutoff <- as.integer(args[8])
Genome1 <- args[9]
libtype <- args[10]

Name1 <- paste(Name, Genome1, sep=".")

File1 <- paste(Name1, ".unfiltered.counts.csv", sep="")
File2 <- paste(Name1, ".filtered.counts.csv", sep="")

# print(paste(dir, Namehg, File1, sep="/"))
hg_unfil <- read.csv(paste(dir, File1, sep="/"), header=T)
hg_fil <- read.csv(paste(dir, File2, sep="/"), header=T)

Df <- cbind(hg_fil, hg_unfil$fragments)

Genome1_unique <- paste(Genome1, "unique", sep ="_")
Genome1_total <- paste(Genome1, "total", sep ="_")
Genome1_dup <- paste(Genome1, "dup", sep ="_")


colnames(Df) <- c("R1", "R2", "R3", "P5", Genome1_unique, Genome1_total)
#head(Df)

#Df$uniqueNuc <- Df$hg_unique + Df$mm_unique
#Df$totalNuc <- Df$hg_total + Df$mm_total
Df$dup <- 1 - Df[ ,5]/Df[ ,6]
colnames(Df)[7] <- Genome1_dup

Genome1_libsize <- paste(Genome1, "libsize", sep ="_")
Df$temp <- 0
colnames(Df)[8] <- Genome1_libsize
#head(Df)

# remove barcode combination with less than N reads
Df <- Df[Df[ ,6] >= Cutoff, ]

for (i in 1:nrow(Df)){
  if (Df[i,6] != Df[i,5]) {
    Df[i,8] <- round(uniroot(function(x) (x*(1-exp(-Df[i,6]/x))-Df[i,5]),
                               lower = 0, upper = 10000000, tol = 1e-7)$root)}
  else{
    Df[i,8] <- Df[i,6]
  }
}
# for pair-end deviced library size by 2 
if (libtype != "RNA"){
   if (libtype != "ATAC"){
      Df[ ,8] <- Df[ ,8]/2
   }
}
#head(Df)
print(paste("Average ", Genome1, " lib size of ", Name, " is:", sep=""))
mean(Df[Df[, 8] > Cutoff,8])

alllibeize <- sum(Df[ ,8])
duprate <- mean(Df[Df[, 8] > 500,7])
allreads <- sum(Df[ ,6])
dupreads <- sum(Df[ ,6]) - sum(Df[ ,5])
out = data.frame(LIBRARY="Unkown",
		 UNPAIRED_READS_EXAMINED=0,
		 READ_PAIRS_EXAMINED=allreads ,
		 SECONDARY_OR_SUPPLEMENTARY_RDS=0,
		 UNMAPPED_READS=0,
		 UNPAIRED_READ_DUPLICATES=0,
		 READ_PAIR_DUPLICATES= dupreads,
		 READ_PAIR_OPTICAL_DUPLICATES=0,
		 PERCENT_DUPLICATION=duprate,
		 ESTIMATED_LIBRARY_SIZE=alllibeize)
fileout <- paste(Name, Genome1, 'dups.log', sep=".")
write.table(out, fileout, quote=F, row.names=F, col.names=T, sep="\t")

# PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
# Brian_ATAC.mm10.dups.log

# Plot lib size
print("Plot Libsize, Species mixing etc...")

Size2 <-Df[, 8];  Size2 <- Size2[Size2 >0]; # plot(Size)
file6 <- paste(Name, Genome1, 'libSize.pdf', sep=".")
pdf(paste(dir,file6, sep="/"))
plot(log10(sort(Size2, decreasing = T)), xlab="Barcode rank", ylab = "log10 (lib size)", main = "Estimated Lib Sizes per Cell", col="darkblue", pch=16)
garbage <- dev.off()

file12 <- paste(Name, Genome1, 'FilvsUnfil.pdf', sep=".")
pdf(paste(dir,file12, sep="/"))
plot(Df[ , 5],Df[, 6], xlab=paste(Genome1, " unique molecules", sep=""), ylab = "total reads", main = "Saturation", col="darkblue", pch=16)
abline(1,1, col="blue")
garbage <- dev.off()

temp <- paste(Name,".counts.csv", sep="")
File <- paste(dir, temp, sep="/")
write.csv(Df, File)

# C/X = 1 - exp( -N/X ) 
# 
# where X = number of distinct molecules in library 
# N = number of read pairs
# C = number of distinct fragments observed in read pairs"

