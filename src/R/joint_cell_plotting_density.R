library(ggplot2)
library(ggpointdensity)

args <- commandArgs()
pkr <- args[6]
barcode_metadata_file <- args[7]

options(scipen=999)

barcode_metadata <- read.csv(barcode_metadata_file)
passing_df <- barcode_metadata[barcode_metadata$QC %in% c("RNA only", "ATAC only", "both"),]

# get max x and y coords to set plot limits
round_to_power_10 <- function(x){
  return(10^ceiling(log10(x)))
}
max_x <- max(passing_df$frags)
max_y <- max(passing_df$umis)
xy_lim <- round_to_power_10(max(max_x, max_y))

# palette from https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
sambaNight <- c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500')

density_plot <- ggplot(passing_df, aes(x=frags, y=umis)) + 
                        geom_pointdensity(size=0.7) +
                        scale_color_gradientn(colors=sambaNight) + 
                        labs(title=paste0("Joint Cell Calling (", pkr, "): Density Plot", sep=""),
                             x="ATAC Unique Fragments per Barcode",
                             y="RNA UMIs per Barcode") +
                        theme_light() +
                        theme(plot.margin=margin(t=9, r=36.5, b=25, l=9, unit="pt"),
                              plot.title=element_text(size=12.5, hjust=0.5),
                              axis.title=element_text(size=11),
                              axis.text=element_text(size=8.5),
                              legend.title=element_text(size=8),
                              legend.text=element_text(size=6),
                              panel.grid.minor=element_blank()) +
                        scale_x_continuous(trans="log10",
                                           limits=c(10,xy_lim)) +
                        scale_y_continuous(trans="log10",
                                           limits=c(10,xy_lim))
                        
out_file <- paste0(pkr, "_joint_cell_density_plot.png", sep="")
# density plotting fails when 0 cells pass both filters
if (sum(barcode_metadata$QC=="both") > 0) {
  png(out_file, width=8.75, height=6, units="in", res=300)
  density_plot
  dev.off()
}
