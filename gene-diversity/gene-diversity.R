library(ggplot2)
library(grid)
library(vegan)
library(corrplot)

#dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/sharing.txt')

in.file <- file.choose()
df <- read.csv(in.file)

in.file <- file.choose()
dfref <- read.csv(in.file)
df <- merge(df,dfref,by="locus")

dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/id-name-cat-conv.csv')
names(dfref)[1] <- "locus"
df <- merge(df,dfref,by="locus")


# setting up variables
df$allelic.div <- c( log10(  df$count.nuc / df$avg.length) )
df$div.diff <- c( (df$count.nuc - (mean(df$count.nuc / df$avg.length) * df$avg.length) ) / df$avg.length )

df$p.vs.nuc <- (df$avg.length - df$varsites.x) / df$avg.length
df$p.vs.aa <- ((df$avg.length / 3) - df$varsites.y) / (df$avg.length / 3)
df$r.count <- df$count.aa / df$count.nuc
df$r.vs <- df$p.vs.aa / df$p.vs.nuc

# removing worst loci in terms of tagging (usually ones with seed problems)
df.all <- df
cutoff <- 5000
df <- df[df$missing < cutoff,]


# correlation plots for all numeric variables
df2 <- subset(df, select = -c(locus,name,category,drug.resistance) )
M <- cor(df2)
corrplot(M, method = "circle")




# For my own peace of mind...
levels(df$category)[levels(df$category)=="Cell wall and cell processes"] <- "Cell processes"
levels(df$category)[levels(df$category)=="Insertion seqs and phages"] <- "IS and phages"
levels(df$category)[levels(df$category)=="Conserved hypotheticals"] <- "Hypotheticals"
levels(df$category)[levels(df$category)=="Lipid metabolism"] <- "Lipids"
levels(df$category)[levels(df$category)=="Information pathways"] <- "Signaling"
levels(df$category)[levels(df$category)=="Intermediary metabolism and respiration"] <- "Metabolism"
levels(df$category)[levels(df$category)=="Regulatory proteins"] <- "Regulatory"
levels(df$category)[levels(df$category)=="Ribosomal protein"] <- "Ribosomal"
levels(df$category)[levels(df$category)=="Virulence, detoxification and adaptation"] <- "Virulence"

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Distribution of missing (vast majority less than 200)
# ggplot(df, aes(x=missing)) + geom_density() + geom_vline(xintercept=200) + xlim(0, 6000)

# that's all the data set up done!


# Total distribution
d1 <- 
  ggplot(df, aes(x=allelic.div, fill="")) + geom_histogram(alpha=.7, binwidth=.025) + 
  ggtitle("Density of alleles per nucleotide position (log transformed)") +
  theme(plot.title = element_text(face="bold"), axis.text.x = element_blank(), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_vline(xintercept=mean(df$allelic.div), colour = "red", linetype="dashed") +
  ylab("Density")  + xlab("Alleles per nucleotide, in log10 scale") + guides(fill=FALSE) 
# Distribution by sharing
m.table <- aggregate(allelic.div ~ sharing2, data = df, mean)      
d2 <- ggplot(df, aes(x=allelic.div, fill=sharing2)) + geom_histogram(alpha=.5, position="identity", binwidth=.025) + 
  ggtitle("Density of alleles per nucleotide position (log transformed), separated by sharing status") +
  theme(plot.title = element_text(face="bold"), axis.text.x = element_blank(), axis.title.x=element_text(vjust=-.5, size=14)) +
  scale_fill_manual(values=c("#56B4E9","#E69F00"), name = "Sharing status") + 
  theme(legend.position=c(.9, .8)) +
  geom_vline(xintercept=m.table[1,2], colour = "blue", linetype="dashed") +
  geom_vline(xintercept=m.table[2,2], colour = "brown", linetype="dashed") +
  ylab("Density")  + xlab("Alleles per nucleotide, in log10 scale") 

multiplot(d1,d2)



p1 <- ggplot(df, aes(x=allelic.div)) + geom_histogram(alpha=.7, binwidth = 0.01)
p2 <- ggplot(df, aes(x=div.diff)) + geom_histogram(alpha=.7, binwidth = 0.001)
p3 <- ggplot(df, aes(x=r.vs)) + geom_histogram(alpha=.7, binwidth = 0.001)
p4 <- ggplot(df, aes(x=r.count)) + geom_histogram(alpha=.7, binwidth = 0.005)

multiplot(p1,p2,p3,p4)


#Scatter plot of diversity by category coloured by sharing, same coloured by missing
ggplot(df) +
  geom_boxplot(aes(x=category, y=allelic.div), alpha = 0, outlier.shape = " ") +
  geom_point( data=subset(df, allelic.div < -.9 & allelic.div > -1.7),aes(x=category, y=allelic.div, colour=drug.resistance), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) + 
  geom_point( data=subset(df, allelic.div > -.9 | allelic.div < -1.7),aes(x=category, y=allelic.div, colour=drug.resistance), size=5, alpha=.5) +
  coord_flip() + scale_fill_manual(values=c("#56B4E9","#E69F00")) + 
  guides(col = guide_legend(title = "Drug Resistance")) +
  theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.text.x = element_blank(), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$allelic.div), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by association with drug resistance") +
  xlab("")  + ylab("Alleles per nucleotide, in log10 scale") +
  geom_text( data=subset(df, allelic.div > -.95 | allelic.div < -1.7), aes(x=category, y=allelic.div, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by length
ggplot(df) +
  geom_boxplot(aes(x=category, y=allelic.div), alpha = 0, outlier.shape = " ") +
  geom_point( data=subset(df, allelic.div < -.9 & allelic.div > -1.7),aes(x=category, y=allelic.div, colour=length), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) + 
  geom_point( data=subset(df, allelic.div > -.9 | allelic.div < -1.7),aes(x=category, y=allelic.div, colour=length), size=5, alpha=.5) +
  coord_flip() +  scale_colour_gradient(limits=c(0,1000), low="red", high="grey") +
  theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.text.x = element_blank(), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$allelic.div), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Alleles per nucleotide, in log10 scale") +
  geom_text( data=subset(df, allelic.div > -.95 | allelic.div < -1.7), aes(x=category, y=allelic.div, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# with null model measure
ggplot(df) +
  geom_boxplot(aes(x=category, y=div.diff), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=div.diff, colour=length), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + scale_colour_gradient(limits=c(0,500), low="red", high="grey") +
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$div.diff), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Difference in alleles/nucleotide from expected based on mean") +
  geom_text( data=subset(df, div.diff > 0.07 | div.diff < -0.03), aes(x=category, y=div.diff, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by Simpson's Diversity Index 
ggplot(df) +
  geom_boxplot(aes(x=category, y=div.index), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=div.index, colour=length), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + scale_colour_gradient(limits=c(0,500), low="red", high="grey") +
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$div.index), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Simpson's diversity index per locus") +
  geom_text( data=subset(df, div.index > 0.07 | div.index < -0.03), aes(x=category, y=div.diff, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by proportion of variable sites
ggplot(df) +
  geom_boxplot(aes(x=category, y=pvarsites.nuc), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=pvarsites.nuc, colour=psynnon), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$pvarsites), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Difference in alleles/nucleotide from expected based on mean") +
  geom_text( data=subset(df, pvarsites.nuc > 0.1 | pvarsites.nuc < 0.005), aes(x=category, y=pvarsites.nuc, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by ratio of unique amino acid sequences to unique nucleotide 
ggplot(df) +
  geom_boxplot(aes(x=category, y=psynnon), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=psynnon, colour=drug.resistance), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$psynnon), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Ratio of unique alleles in nucleotide to amino acid format") +
  geom_text( data=subset(df, psynnon > .8 | psynnon < .4), aes(x=category, y=psynnon, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 


ggplot(df) +
  geom_boxplot(aes(x=category, y=ratio.varsites), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=ratio.varsites, colour=drug.resistance), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$ratio.varsites), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Ratio of unique alleles in nucleotide to amino acid format") +
  geom_text( data=subset(df, ratio.varsites > .8 | ratio.varsites < .4), aes(x=category, y=ratio.varsites, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 




ggplot(df, aes(pvarsites, allelic.div)) + geom_point() +
  geom_text( data=subset(df, pvarsites < 0.01), aes(x=pvarsites, y=allelic.div, label=name), size=4)

ggplot(df, aes(pvarsites, div.diff)) + geom_point() +
  geom_text( data=subset(df, div.diff > 0.15), aes(x=pvarsites, y=div.diff, label=name), size=4)

