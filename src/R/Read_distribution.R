#!/usr/bin/Rscript

args <- commandArgs()
#print(args)

dir <- args[6]
Name <- args[7]

#tail -n +5 RNA_2500.mm10.read_distribution.txt | head -n -1 > temp2.txt
#head -n 3 RNA_2500.mm10.read_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g'

Df <- read.table(paste(dir,"temp1.txt",sep="/"), header = T)

Total_tag <- as.numeric(read.table(paste(dir,"temp2.txt",sep="/"), header = F)); Total_tag
Df2 <- data.frame(Df$Group,Df$Tag_count)
colnames(Df2) <- c("Group", "Tag_count")
Inter <- data.frame("Intergenic",Total_tag-sum(Df2$Tag_count[1:4]))
colnames(Inter) <- c("Group", "Tag_count")
Df3 <- rbind(Df2,Inter)
Df3$Tag_perc <- Df3$Tag_count/Total_tag*100; Df3

temp <- paste(Name,'.reads_distribution.pdf', sep="")

pdf(paste(dir,temp, sep="/"))
barplot(Df3$Tag_perc[c(1:4,11)], names.arg=Df3$Group[c(1:4,11)], cex.names = 0.9,
        ylab="Percentage of reads (%)", main="Reads distribution")
dev.off()

temp <- paste(Name,".reads_distribution2.txt", sep="")
write.table(Df3,paste(dir, temp, sep="/"), quote = F, row.names=F, sep="\t")