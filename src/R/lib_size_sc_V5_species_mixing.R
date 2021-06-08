#!/usr/bin/Rscript

args <- commandArgs()
#print(args)

dir <- args[6]
Name <- args[7]
Cutoff <- as.integer(args[8])
libtype <- args[9]

Namehg <- paste(Name, "hg19", sep=".")
Namemm <- paste(Name, "mm10", sep=".")

File1 <- paste(Namehg, ".unfiltered.counts.csv", sep="")
File2 <- paste(Namehg, ".filtered.counts.csv", sep="")
File3 <- paste(Namemm, ".unfiltered.counts.csv", sep="")
File4 <- paste(Namemm, ".filtered.counts.csv", sep="")

# print(paste(dir, Namehg, File1, sep="/"))
hg_unfil <- read.csv(paste(dir, File1, sep="/"), header=T)
hg_fil <- read.csv(paste(dir, File2, sep="/"), header=T)
mm_unfil <- read.csv(paste(dir, File3, sep="/"), header=T)
mm_fil <- read.csv(paste(dir, File4, sep="/"), header=T)


Df <- cbind(hg_fil,mm_fil$fragments, hg_unfil$fragments, mm_unfil$fragments)
colnames(Df) <- c("R1", "R2", "R3", "P5", "hg_unique",  "mm_unique", "hg_total", "mm_total")
#head(Df)
Df$uniqueNuc <- Df$hg_unique + Df$mm_unique
Df$totalNuc <- Df$hg_total + Df$mm_total
Df$dup <- 1 - Df$uniqueNuc/Df$totalNuc 
#head(Df)
#dim(Df)
# remove barcode combination with less than N reads
Df <- Df[Df$totalNuc >= Cutoff, ]
Df$libsize <- 0
# dim(Df)
# head(Df)
for (i in 1:nrow(Df)){
  if (Df$totalNuc[i] == 0 ) {
     Df$libsize[i] <- 0
   } else if (Df$totalNuc[i] != Df$uniqueNuc[i]) {
    Df$libsize[i] <- round(uniroot(function(x) (x*(1-exp(-Df$totalNuc[i]/x))-Df$uniqueNuc[i]), 
                               lower = 0, upper = 10000000, tol = 1e-7)$root)}
  else{
    Df$libsize[i] <- Df$totalNuc[i]
  }
}
if (libtype != "RNA"){
   if (libtype != "ATAC"){
      Df$libsize <- Df$libsize/2
   }
}
head(Df)
print(paste("Average lib size of ", Name, " is:", sep=""))
mean(Df[Df$libsize > Cutoff,]$libsize)

Df$hg_libsize <- 0
for (i in 1:nrow(Df)){
  if (Df$hg_total[i] == 0 ) {
     Df$hg_libsize[i] <- 0
  } else if (Df$hg_total[i] != Df$hg_unique[i]) {
    Df$hg_libsize[i] <- round(uniroot(function(x) (x*(1-exp(-Df$hg_total[i]/x))-Df$hg_unique[i]),
                               lower = 0, upper = 10000000, tol = 1e-7)$root)}
  else{
    Df$hg_libsize[i] <- Df$hg_total[i]
  }
}
if (libtype != "RNA"){
   Df$hg_libsize <- Df$hg_libsize/2
}

# head(Df)

print(paste("Average hg19 lib size of ", Name, " is :", sep=""))
mean(Df[Df$hg_libsize > 500,]$hg_libsize)

Df$mm_libsize <- 0
for (i in 1:nrow(Df)){
  if (Df$mm_total[i] == 0 ) {
     Df$mm_libsize[i] <- 0
     } else if (Df$mm_total[i] != Df$mm_unique[i]) {
    Df$mm_libsize[i] <- round(uniroot(function(x) (x*(1-exp(-Df$mm_total[i]/x))-Df$mm_unique[i]),
                               lower = 0, upper = 10000000, tol = 1e-7)$root)}
  else{
    Df$mm_libsize[i] <- Df$mm_total[i]
  }
}
if (libtype != "RNA"){
   Df$mm_libsize <- Df$mm_libsize/2
}

# head(Df)
print(paste("Average mm10 lib size of ", Name, " is :", sep=""))
mean(Df[Df$mm_libsize > 500,]$mm_libsize)

# make a dups.log file
alllibeize <- sum(Df$hg_libsize)
duprate <- mean(Df[Df$hg_libsize > 500,]$dup)
allreads <- sum(Df$hg_total)
dupreads <- sum(Df$hg_total) - sum(Df$hg_unique)
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
fileout <- paste(Name, 'hg19.dups.log', sep=".")
write.table(out, fileout, quote=F, row.names=F, col.names=T)

alllibeize <- sum(Df$mm_libsize)
duprate <- mean(Df[Df$mm_libsize > 500,]$dup)
allreads <- sum(Df$mm_total)
dupreads <- sum(Df$mm_total) - sum(Df$mm_unique)
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
fileout <- paste(Name, 'mm10.dups.log', sep=".")
write.table(out, fileout, quote=F, row.names=F, col.names=T)

# Plot lib size
#Size <-Df$libsize;  Size <- Size[Size >0]; # plot(Size)
#file5 <- paste(Name,'.libSize.pdf', sep="")
#pdf(paste(dir,file5, sep="/"))
#plot(log10(sort(Size, decreasing = T)), xlab="Barcode rank", ylab = "log10 (lib size)", main = "Estimated Lib Sizes per Cell", col="darkred", pch=16)
#garbage <- dev.off()

print("Plot Libsize, Species mixing etc...")

Size2 <-Df$hg_libsize;  Size2 <- Size2[Size2 >0]; # plot(Size)
file6 <- paste(Name,'.hg19.libSize.pdf', sep="")
pdf(paste(dir,file6, sep="/"))
plot(log10(sort(Size2, decreasing = T)), xlab="Barcode rank", ylab = "log10 (lib size)", main = "Estimated Lib Sizes per Cell", col="darkblue", pch=16)
garbage <- dev.off()

Size3 <-Df$mm_libsize;  Size3 <- Size3[Size3 >0]; # plot(Size)
file7 <- paste(Name,'.mm10.libSize.pdf', sep="")
pdf(paste(dir,file7, sep="/"))
plot(log10(sort(Size3, decreasing = T)), xlab="Barcode rank", ylab = "log10 (lib size)", main = "Estimated Lib Sizes per Cell", col="orange", pch=16)
garbage <- dev.off()

file8 <- paste(Name,'.mixing.log10.libSize.pdf', sep="")
pdf(paste(dir,file8, sep="/"))
plot(log10(Df$hg_libsize+1),log10(Df$mm_libsize+1), xlab="hg19 log10 (lib size)", ylab = "mm10 log10 (lib size)", main = "Estimated Lib Sizes per Cell", col="forestgreen", pch=16)
garbage <- dev.off()

file9 <- paste(Name,'.mixing.normal.libSize.pdf', sep="")
pdf(paste(dir,file9, sep="/"))
plot(Df$hg_libsize,Df$mm_libsize, xlab="hg19 lib size", ylab = "mm10 lib size", main = "Estimated Lib Sizes per Cell", col="forestgreen", pch=16)
garbage <- dev.off()

file10 <- paste(Name,'.mixing.log10.uniFrag.pdf', sep="")
pdf(paste(dir,file10, sep="/"))
plot(log10(Df$hg_unique+1),log10(Df$mm_unique+1), xlab="hg19 log10 (unique molecules)", ylab = "mm10 log10 (unique molecules)", main = "Unique molecules per Cell", col="darkred", pch=16)
garbage <- dev.off()

file11 <- paste(Name,'.mixing.normal.uniFrag.pdf', sep="")
pdf(paste(dir,file11, sep="/"))
plot(Df$hg_unique,Df$mm_unique, xlab="hg19 unique molecules", ylab = "mm10 unique molecules", main = "Unique molecules per Cell", col="darkred", pch=16)
garbage <- dev.off()

file12 <- paste(Name,'.hg19.FilvsUnfil.pdf', sep="")
pdf(paste(dir,file12, sep="/"))
plot(Df$hg_unique,Df$hg_total, xlab="hg19 unique molecules", ylab = "hg19 total reads", main = "Saturation", col="darkblue", pch=16)
abline(1,1, col="blue")
garbage <- dev.off()

file13 <- paste(Name,'.mm10.FilvsUnfil.pdf', sep="")
pdf(paste(dir,file13, sep="/"))
plot(Df$mm_unique,Df$mm_total, xlab="mm10 unique molecules", ylab = "mm10 total reads", main = "Saturation", col="orange", pch=16)
garbage <- dev.off()

temp <- paste(Name,".counts.csv", sep="")
File <- paste(dir, temp, sep="/")
write.csv(Df, File)

# C/X = 1 - exp( -N/X ) 
# 
# where X = number of distinct molecules in library 
# N = number of read pairs
# C = number of distinct fragments observed in read pairs"

