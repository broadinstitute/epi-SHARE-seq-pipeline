#!/usr/bin/Rscript

args <- commandArgs()
#print(args)

dir <- args[6]
Name <- args[7]

# read csv
Counts <- read.csv(paste(dir, Name, sep="/"), header=F, sep = "")
colnames(Counts) <- c("Freq","Barcode")
#head(Table)

#split barcode name
Counts$P5 <- substr(Counts$Barcode,19,23)
#head(Counts)
P5s <- sort(unique(Counts$P5))
print("P5s are"); print(P5s)

# generate all combinations
R1s <- paste("R1.", formatC(1:96, width=2, flag="0"), sep="")
R2s <- paste("R2.", formatC(1:96, width=2, flag="0"), sep="")
R3s <- paste("R3.", formatC(1:96, width=2, flag="0"), sep="")

library(tidyr)
Df <- crossing(P5 = P5s, R1 = R1s, R2 = R2s, R3 =  R3s)
Df$Barcode  <- paste(Df$R1, Df$R2, Df$R3, Df$P5, sep=",")
#Df

# match 
Df$Count <- Counts$Freq[match(Df$Barcode, Counts$Barcode)]
Df[is.na(Df)] <- 0
head(Df[order(-Df$Count),])

#stop("")

# re-organize format
library(tibble)
Df2 <- tibble(R1=Df$R1, R2=Df$R2, R3=Df$R3, P5=Df$P5, fragments=Df$Count)
#head(Df2)

temp <- paste(Name, ".csv", sep="")
File <- paste(dir, temp, sep="/")
write.csv(Df2, File, quote = F, row.names = F)
print(paste("Finished counting for ", Name, sep=""))

#stop()
