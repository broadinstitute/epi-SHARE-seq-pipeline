#!/usr/bin/Rscript

args <- commandArgs()
#print(args)

inp <- args[6]
out <- args[7]

# read csv
Counts <- read.csv(inp, header=F, sep = "")
colnames(Counts) <- c("Freq","Barcode")
#head(Table)

#split barcode name
#Counts$P5 <- substr(Counts$Barcode,19,23)
#head(Counts)
#P5s <- sort(unique(Counts$P5))
#print("P5s are"); print(P5s)


barcodes=unlist(lapply(Counts$Barcode,function(x){strsplit(x,",")[[1]][c(1,2,3)]}))
print(length(unique(barcodes)))
#If empty is giving NA
pkrs=unlist(lapply(Counts$Barcode,function(x){strsplit(x,",")[[1]][4]}))
# generate all combinations
#R1s <- paste("R1.", formatC(1:96, width=2, flag="0"), sep="")
#R2s <- paste("R2.", formatC(1:96, width=2, flag="0"), sep="")
#R3s <- paste("R3.", formatC(1:96, width=2, flag="0"), sep="")

library(tidyr)
library(dplyr)

Df <- crossing(R1 = barcodes, R2 = barcodes, R3 = barcodes, PKR = pkrs)
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
