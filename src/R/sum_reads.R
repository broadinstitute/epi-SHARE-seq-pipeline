#!/usr/bin/Rscript
# Mei will fix this
args <- commandArgs()
#print(args)

inp <- args[6]
out <- args[7]
barcodes_list_fnp <- args[8]

# Read barcodes
observed_barcodes = read.csv(barcodes_list_fnp, header=F, sep = "")
# Convert into one line
barcodes.R1=unlist(lapply(observed_barcodes$V1,function(x){strsplit(x,",")[[1]][c(1)]}))
barcodes.R2=unlist(lapply(observed_barcodes$V1,function(x){strsplit(x,",")[[1]][c(2)]}))
barcodes.R3=unlist(lapply(observed_barcodes$V1,function(x){strsplit(x,",")[[1]][c(3)]}))

# read csv
Counts <- read.csv(inp, header=F, sep = "")
colnames(Counts) <- c("Freq","Barcode")
#head(Table)

#split barcode name
#Counts$P5 <- substr(Counts$Barcode,19,23)
#head(Counts)
#P5s <- sort(unique(Counts$P5))
#print("P5s are"); print(P5s)


#barcodes=unlist(lapply(Counts$Barcode,function(x){strsplit(x,",")[[1]][c(1,2,3)]}))
print(length(unique(barcodes.R1)))
print(length(unique(barcodes.R2)))
print(length(unique(barcodes.R3)))
#If empty is giving NA
pkrs=unlist(lapply(Counts$Barcode,function(x){strsplit(x,",")[[1]][4]}))
# generate all combinations
#R1s <- paste("R1.", formatC(1:96, width=2, flag="0"), sep="")
#R2s <- paste("R2.", formatC(1:96, width=2, flag="0"), sep="")
#R3s <- paste("R3.", formatC(1:96, width=2, flag="0"), sep="")

library(tidyr)
library(dplyr)

Df <- crossing(R1 = barcodes.R1, R2 = barcodes.R2, R3 = barcodes.R3, PKR = pkrs)
if(is.na(Df$PKR[1])){
    Df = dplyr::select(Df,-PKR)
}
Df$Barcode <- apply(Df, 1, paste, collapse = ",")

# match
Df$fragments <- Counts$Freq[match(Df$Barcode, Counts$Barcode)]
Df[is.na(Df)] <- 0
head(Df[order(-Df$fragments),])

Df=dplyr::select(Df,-Barcode)

#stop("")

# re-organize format
#library(tibble)
#Df2 <- tibble(R1=Df$R1, R2=Df$R2, R3=Df$R3, fragments=Df$Count)
#head(Df2)

write.csv(Df, out, quote = F, row.names = F)

#stop()
