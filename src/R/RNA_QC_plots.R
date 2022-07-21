## Make pngs of detected genes vs barcode rank plot, detected genes vs UMIs plot,
## and UMIs vs barcode rank plot from h5 matrix

args <- commandArgs(trailingOnly=T)
input_file <- args[1] #h5 RNA count matrix 

suppressMessages(library(Seurat))

# import h5 matrix
h5_matrix <- Read10X_h5(input_file)

print(ncol(h5_matrix))

UMIs_per_barcode <- colSums(h5_matrix)
genes_per_barcode <- colSums(h5_matrix>0)

print("Top 10 UMI counts: ")
head(sort(UMIs_per_barcode, decreasing=T), 10)

df <- data.frame(UMIs_per_barcode, genes_per_barcode)
# sort df in descending order of genes_per_barcode
df <- df[order(-genes_per_barcode),]
# add column containing barcode rank (based on number of genes)
df$barcode_rank <- c(1:nrow(df))

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

# plot log10(detected genes) vs barcode rank
png("out.genes.vs.barcode.rank.png")
elbow <- elbow_finder(df$barcode_rank, log10(df$genes_per_barcode))
plot(x=df$barcode_rank, y=log10(df$genes_per_barcode),
     xlab="Barcode Rank", ylab="log10(Detected Genes)",
     main="Detected Genes vs Barcode Rank",
     col="darkblue",
     pch=16)
abline(v=elbow[1], h=elbow[2])
text(elbow[1], elbow[2],
     paste(elbow[1], round(elbow[2], digits=2), sep=", "),
     adj=c(-0.1,-0.5))
legend("topright", c(paste("Number of barcodes:              ", nrow(df), sep=""),
                    paste("Barcodes with > 10 genes:     ", sum(df$genes_per_barcode>10), sep=""),
                    paste("Barcodes with > 100 genes:   ", sum(df$genes_per_barcode>100), sep=""),
                    paste("Barcodes with > 1000 genes: ", sum(df$genes_per_barcode>1000), sep ="")),
       bty="n")
dev.off()

# plot detected genes vs UMIs
png("out.genes.vs.UMIs.png")
plot(x=df$UMIs_per_barcode, y=df$genes_per_barcode,
     xlab="UMIs", ylab="Detected Genes",
     main="Detected Genes vs UMIs",
     col="darkblue",
     pch=16)
dev.off()

# plot UMIs vs barcode rank
png("out.UMIs.vs.barcode.rank.png")
elbow <- elbow_finder(df$barcode_rank, df$UMIs_per_barcode)
plot(x=df$barcode_rank, y=df$UMIs_per_barcode,
     xlab="Barcode Rank", ylab="UMIs",
     main="UMIs vs Barcode Rank",
     col="darkblue",
     pch=16)
abline(v=elbow[1], h=elbow[2])
text(elbow[1], elbow[2],
     paste(elbow[1], elbow[2], sep=", "),
     adj=c(-0.1,-0.5))
legend("topright", c(paste("Number of barcodes:               ", nrow(df), sep=""),
                     paste("Barcodes with > 10 UMIs:       ", sum(df$UMIs_per_barcode>10), sep=""),
                     paste("Barcodes with > 100 UMIs:     ", sum(df$UMIs_per_barcode>100), sep=""),
                     paste("Barcodes with > 1000 UMIs:   ", sum(df$UMIs_per_barcode>1000), sep="")),
       bty="n")
dev.off()
