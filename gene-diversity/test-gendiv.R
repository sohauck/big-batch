Args <- commandArgs(TRUE);      # retrieve args

print(Args[1])

# Setting up packages
library(ggplot2)

# Going to directory with results
setwd("/Volumes/sofia/big-batch/gene-diversity/examples/GeneDiv-MYCO-table-F/")

# Reading in table with output from Perl runs
df <- read.table("ResultsTable.txt", header = TRUE, sep = "\t")

# Setting up new variables 
df$AllelicDiv <- c( df$CountNuc / df$AvgLength )
df$ADivNM <- c( (df$CountNuc - (mean(df$CountNuc / df$AvgLength) * df$AvgLength) ) / df$AvgLength )

df$VSitesNuc <- (df$AvgLength - df$NonVarNuc) / df$AvgLength
df$VSitesAA <- ((df$AvgLength / 3) - df$NonVarAA) / (df$AvgLength / 3)

df$RatioCount <- df$CountAA / df$CountNuc
df$RatioVS <- df$VSitesAA / df$VSitesNuc

#df$allelic.div.z <- scale(df$allelic.div, center = TRUE, scale = TRUE)


# Writing fuller table back onto results directory
write.table (df, "CalculationsTable.txt", sep = "\t", row.names = FALSE)

# Moving to Graphs folder for saving graphs
setwd("/Volumes/sofia/big-batch/gene-diversity/examples/GeneDiv-MYCO-table-F/Graphs")

# Making graphs!

# Histogram of loci with >90% missing loci (actually all, but with dashed line)
ggplot(df, aes(x=Missing)) + 
  geom_histogram(binwidth=(nrow(df)/30)) +
  geom_vline(xintercept=(9/10)*nrow(df),size=1, colour="red", linetype="dashed") + 
  theme_minimal() + scale_x_continuous(limits=c(0, nrow(df)),breaks=c((1:4*(1/4))*nrow(df)),labels=c("75%","50%","25%","00%")) +
  theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=14), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  ggtitle("Distribution of loci by percentage of isolates where they have a known sequence") +
  ylab("") + xlab("Percentage of isolates in which locus is known")

ggsave("Dist-Missing.png", height = 9, width = 12, dpi = 100)

# Distribution of AllelicDiv and ADivNM
ggplot(df, aes(x=AllelicDiv, fill="")) + 
  geom_histogram(alpha=.7) + 
  ggtitle("Density of alleles per nucleotide position") +
  theme(plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_vline(xintercept=mean(df$AllelicDiv), colour = "red", linetype="dashed") +
  ylab("Density")  + xlab("Alleles per length in nucleotides") + guides(fill=FALSE) 

ggsave("Dist-AllelicDiv.png", height = 9, width = 12, dpi = 100)


ggplot(df, aes(x=ADivNM, fill="")) + 
  geom_histogram(alpha=.7) + 
  ggtitle("Density of alleles per nucleotide position, compared to null model") +
  theme(plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_vline(xintercept=mean(df$ADivNM), colour = "red", linetype="dashed") +
  ylab("Density")  + xlab("Difference in alleles per length in nucleotides from genome average") + guides(fill=FALSE) 

ggsave("Dist-ADivNM.png", height = 9, width = 12, dpi = 100)


# The famous bubble plot

dodge3 <- position_dodge(width = 0.5)

set.seed(20061001, "Mersenne-Twister")

set.seed(123)
x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE)

ggplot(df) +
  geom_point( data=df , aes(x=x1, y=ADivNM), size=5, alpha=.5, label=df$Locus) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_blank(), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14), axis.text.y = element_blank()) +
  geom_hline(yintercept=mean(df$ADivNM), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci") + xlab("")  + ylab("Alleles per nucleotide, by z-score") +
  geom_text( data=df, aes(x=x1, y=ADivNM, label=Locus), size=4, alpha=.8, vjust=-1.5) 





ggplot(df) +
  geom_point( data=subset(df, allelic.div.z < 1 | allelic.div.z > -1),aes(x=1, y=allelic.div.z), size=5, alpha=.5, position = position_jitter(w=0.5, h=0)) + 
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_blank(), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14), axis.text.y = element_blank()) +
  geom_hline(yintercept=mean(df$allelic.div.z), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci") +
  xlab("")  + ylab("Alleles per nucleotide, by z-score") +
  geom_text( data=subset(df, allelic.div.z > 6), aes(x=1, y=allelic.div.z, label=locus), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.5)) 





ggsave("Point-AllelicDiv.png", height = 9, width = 12, dpi = 100)






