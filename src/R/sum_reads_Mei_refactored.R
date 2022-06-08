#!/usr/bin/Rscript
# Generates matrix containing read frequency of all barcode combinations

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

args <- commandArgs()

read_freq_file <- args[6] # read frequency file (frequency, barcodes)
output_file <- args[7] # output filename (R1, R2, R3, PKR, fragments)
observed_barcodes_file <- args[8] # barcode combination file (sorted barcodes w no dups)

# Read in barcodes file
observed_barcodes <- read.csv(observed_barcodes_file, header=F, sep="")
# Get R1, R2, R3 vectors
R1_barcodes <- unlist(lapply(observed_barcodes$V1,function(x){strsplit(x,",")[[1]][1]}))
R2_barcodes <- unlist(lapply(observed_barcodes$V1,function(x){strsplit(x,",")[[1]][2]}))
R3_barcodes <- unlist(lapply(observed_barcodes$V1,function(x){strsplit(x,",")[[1]][3]}))

# Read in read frequency file
read_freq <- read.csv(read_freq_file, header=F, sep="")
colnames(read_freq) <- c("freq", "barcode")
# Get PKR vector
PKRs <- unlist(lapply(read_freq$barcode, function(x){strsplit(x,",")[[1]][4]}))

# Generate dataframe containing all combinations of R1, R2, R3, and PKR
combos <- tidyr::crossing(R1=R1_barcodes, R2=R2_barcodes, R3=R3_barcodes, PKR=PKRs)
# Drop PKR column if NA
if(is.na(combos$PKR[1])){
  combos <- dplyr::select(combos,-PKR)
}
# Concatenate elements to get barcode
combos$barcode <- apply(combos, 1, paste, collapse=",")

# Match possible combinations to observed combinations and get read frequency of matches
combos$fragments <- read_freq$freq[match(combos$barcode, read_freq$barcode)]
combos[is.na(combos)] <- 0
combos <- dplyr::select(combos,-barcode)

# Write output file
write.csv(combos, output_file, quote=F, row.names=F)