Args <- commandArgs(TRUE);      # retrieve args

# Setting up packages
library(ggplot2)

# Going to directory with results
setwd(Args[1])

# Reading in table with output from Perl runs
df <- read.table("ResultsTable.txt", header = TRUE, sep = "\t")

# Setting up new variables 
df$AllelicDiv <- c( df$CountNuc / df$AvgLength )
df$ADivNM <- c( (df$CountNuc - (mean(df$CountNuc / df$AvgLength) * df$AvgLength) ) / df$AvgLength )

df$VSitesNuc <- (df$AvgLength     - df$NonVarNuc) /  df$AvgLength
df$VSitesAA <- ((df$AvgLength / 3) - df$NonVarAA) / (df$AvgLength / 3)

df$RatioCount <- df$CountAA  /  df$CountNuc
df$RatioVS    <- df$VSitesAA /  df$VSitesNuc

#df$allelic.div.z <- scale(df$allelic.div, center = TRUE, scale = TRUE)


# Writing fuller table back onto results directory
write.table (df, "CalculationsTable.txt", sep = "\t", row.names = FALSE)

# Moving to Graphs folder for saving graphs
setwd( paste( Args[1] , "/Graphs/", sep = "" ) )

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
  geom_histogram(alpha=.7, binwidth=((1/30)*(range(df$AllelicDiv)[2]-range(df$AllelicDiv)[1]))) + 
  ggtitle("Density of alleles per nucleotide position") + 
  theme_minimal() + 
  theme(plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_vline(xintercept=mean(df$AllelicDiv), colour = "red", linetype="dashed") +
  ylab("Density")  + xlab("Alleles per length in nucleotides") + guides(fill=FALSE) 

ggsave("Dist-AllelicDiv.png", height = 9, width = 12, dpi = 100)


ggplot(df, aes(x=ADivNM, fill="")) + 
  geom_histogram(alpha=.7, binwidth=((1/30)*(range(df$AllelicDiv)[2]-range(df$AllelicDiv)[1]))) + 
  ggtitle("Density of alleles per nucleotide position, compared to null model") +
  theme_minimal() + 
  theme(plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_vline(xintercept=mean(df$ADivNM), colour = "red", linetype="dashed") +
  ylab("Density")  + xlab("Difference in alleles per length in nucleotides from genome average") + guides(fill=FALSE) 

ggsave("Dist-ADivNM.png", height = 9, width = 12, dpi = 100)


# The famous bubble plots!

x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE)


ggplot(df) +
  geom_point( data=df , aes(x=x1, y=AllelicDiv), size=5, alpha=.5, label=df$Locus) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_blank(), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14), axis.text.y = element_blank()) +
  geom_hline(yintercept=mean(df$AllelicDiv), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci") + xlab("")  + ylab("Alleles per nucleotide") +
  geom_text(aes(x=x1, y=AllelicDiv, label=ifelse(AllelicDiv<head(sort(df$AllelicDiv),10)[10] | AllelicDiv>tail(sort(df$AllelicDiv),10)[1],
                                          as.character(Locus),'')), size=4, alpha=.8, vjust=-1.5) 

ggsave("Point-AllelicDiv.png", height = 9, width = 12, dpi = 100)

x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE)

ggplot(df) +
  geom_point( data=df , aes(x=x1, y=ADivNM), size=5, alpha=.5, label=df$Locus) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_blank(), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14), axis.text.y = element_blank()) +
  geom_hline(yintercept=mean(df$ADivNM), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, compared to null model") + xlab("")  + ylab("Difference in alleles per nucleotide from genome average") +
  geom_text(aes(x=x1, y=ADivNM, label=ifelse(ADivNM<head(sort(df$ADivNM),10)[10] | ADivNM>tail(sort(df$ADivNM),10)[1],
                                             as.character(Locus),'')), size=4, alpha=.8, vjust=-1.5) 

ggsave("Point-ADivNM.png", height = 9, width = 12, dpi = 100)


# Ratio of Count
x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE)

ggplot(df) +
  geom_point( data=df , aes(x=x1, y=RatioCount), size=5, alpha=.5, label=df$Locus) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_blank(), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14), axis.text.y = element_blank()) +
  geom_hline(yintercept=mean(df$RatioCount), size=1, colour="blue", linetype="dashed") +
  ggtitle("Distribution of purifying to diversifying selection") + xlab("")  + ylab("Ratio of unique nucleotide to unique amino acid sequences per locus") +
  geom_text(aes(x=x1, y=RatioCount, label=ifelse(RatioCount<head(sort(df$RatioCount),10)[10] | RatioCount>=tail(sort(df$RatioCount),10)[1],
                                             as.character(Locus),'')), size=4, alpha=.8, vjust=-1.5) 

ggsave("Point-RatioCount.png", height = 9, width = 12, dpi = 100)


# Ratio of variable sites
x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE)

ggplot(df) +
  geom_point( data=df , aes(x=x1, y=RatioVS), size=5, alpha=.5, label=df$Locus) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_blank(), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14), axis.text.y = element_blank()) +
  geom_hline(yintercept=mean(df$RatioVS), size=1, colour="blue", linetype="dashed") +
  ggtitle("Distribution of purifying to diversifying selection") + xlab("")  + ylab("Ratio of variable sites in nucleotide to amino acid format per locus") +
  geom_text(aes(x=x1, y=RatioVS, label=ifelse(RatioVS<head(sort(df$RatioVS),10)[10] | RatioVS>=tail(sort(df$RatioVS),10)[1],
                                                 as.character(Locus),'')), size=4, alpha=.8, vjust=-1.5) 

ggsave("Point-RatioVS.png", height = 9, width = 12, dpi = 100)
