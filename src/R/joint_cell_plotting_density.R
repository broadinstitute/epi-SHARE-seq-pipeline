library(ggplot2)
library(ggpointdensity)

args <- commandArgs()
pkr <- args[6]
barcode_metadata_file <- args[7]
min_UMIs <- args[8]
min_genes <- args[9]
min_TSS <- args[10]
min_frags <- args[11]

options(scipen=999)

barcode_metadata <- read.csv(barcode_metadata_file)
passing_df <- barcode_metadata[barcode_metadata$QC %in% c("RNA only", "ATAC only", "both"),]

# get max x and y coords to set plot limits
round_to_power_10 <- function(x){
  return(10^ceiling(log10(x)))
}
max_x <- max(passing_df$frags)
max_y <- max(passing_df$UMIs)
xy_lim <- round_to_power_10(max(max_x, max_y))

# palette from https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
sambaNight <- c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500')

density_plot <- ggplot(passing_df, aes(x=frags, y=UMIs)) + 
                        geom_pointdensity() +
                        scale_color_gradientn(colors=sambaNight) + 
                        labs(title=paste0("Joint Cell Calling (", pkr, "): Density Plot", sep=""),
                             x="ATAC Unique Fragments per Barcode",
                             y="RNA UMIs per Barcode",
                             caption=paste0("ATAC cutoffs: TSS ≥ ", min_TSS, ", frags ≥ ", min_frags, ". RNA cutoffs: UMIs ≥ ", min_UMIs, ", genes ≥ ", min_genes)) +
                        theme_light() +
                        theme(plot.title=element_text(hjust=0.5),
                              plot.caption=element_text(hjust=0.5),
                              panel.grid.minor=element_blank()) +
                        scale_x_continuous(trans="log10",
                                           limits=c(10,xy_lim)) +
                        scale_y_continuous(trans="log10",
                                           limits=c(10,xy_lim))
                        
out_file <- paste0(pkr, "_joint_cell_density_plot.png", sep="")
png(out_file, width=8, height=6, units="in", res=300)
density_plot
dev.off()
