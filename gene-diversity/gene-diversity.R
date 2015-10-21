library(ggplot2)
library(grid)
library(corrplot)

# reading in data
df <- read.csv(file.choose())

dfref <- read.csv(file.choose())
df <- merge(df,dfref,by="locus")

dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/id-name-cat-conv.csv')
names(dfref)[1] <- "locus"
df <- merge(df,dfref,by="locus")

#dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/sharing.txt')


# setting up variables
df$allelic.div.old <- c( log10(  df$count.nuc / df$avg.length) )
df$allelic.div <- c( df$count.nuc / df$avg.length )
df$div.diff <- c( (df$count.nuc - (mean(df$count.nuc / df$avg.length) * df$avg.length) ) / df$avg.length )

df$p.vs.nuc <- (df$avg.length - df$varsites.x) / df$avg.length
df$p.vs.aa <- ((df$avg.length / 3) - df$varsites.y) / (df$avg.length / 3)

df$r.count <- df$count.aa / df$count.nuc
df$r.vs <- df$p.vs.aa / df$p.vs.nuc


# removing worst loci in terms of tagging (usually ones with seed problems)
df <- df.all
cutoff <- (9/10)*nrow(df.all)
df <- df[df$missing < cutoff,]
# Distribution of missing (vast majority less than 200)
ggplot(df.all, aes(x=missing)) + 
  geom_histogram(binwidth=20) +
  geom_vline(xintercept=(1/10)*nrow(df.all),size=1, colour="red", linetype="dashed") + 
  theme_minimal() + scale_x_continuous(limits=c(0, nrow(df.all)),breaks=c((1:4*(1/4))*nrow(df.all)),labels=c("75%","50%","25%","00%")) +
  theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=14), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  ggtitle("Distribution of number of isolates in which loci are tagged") +
  ylab("") + xlab("Number of isolates in which locus is not tagged")



# correlation plots for all numeric variables
df.num <- subset(df, select = -c(locus,name,category,drug.resistance) )
M <- cor(df.num)
corrplot(M, method = "circle")

# correlation before missing are removed (see how missing correlated negatively with diversity in this case)
df.num0 <- subset(df.all, select = -c(locus,name,category,drug.resistance) )
M0 <- cor(df.num0)
corrplot(M0, method = "circle")



df.mea <- subset(df, select = c(allelic.div,div.diff,p.vs.nuc,p.vs.aa,r.count,r.vs) )
M2 <- cor(df.mea)
corrplot(M2, method = "circle")

# just to keep labels shorter
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



# that's all the data set up done!


# Total distribution
d1 <- 
  ggplot(df, aes(x=allelic.div, fill="")) + geom_histogram(alpha=.7, binwidth=.01) + 
  ggtitle("Density of alleles per nucleotide position") +
  theme(plot.title = element_text(face="bold"), axis.text.x = element_blank(), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_vline(xintercept=mean(df$allelic.div), colour = "red", linetype="dashed") +
  ylab("Density")  + xlab("Alleles per length in nucleotides") + guides(fill=FALSE) 



p1 <- ggplot(df, aes(x=allelic.div2)) + geom_histogram(alpha=.7, binwidth = 0.005)
p2 <- ggplot(df, aes(x=div.diff)) + geom_histogram(alpha=.7, binwidth = 0.003)
p3 <- ggplot(df, aes(x=p.vs.nuc)) + geom_histogram(alpha=.7, binwidth = 0.003)
p4 <- ggplot(df, aes(x=p.vs.aa)) + geom_histogram(alpha=.7, binwidth = 0.005)
p5 <- ggplot(df, aes(x=r.count)) + geom_histogram(alpha=.7, binwidth = 0.004)
p6 <- ggplot(df, aes(x=r.vs)) + geom_histogram(alpha=.7, binwidth = 0.05)

multiplot(p1,p2,p3,p4,p5,p6)


#Scatter plot of diversity by category coloured by drug resistance, 
ggplot(df) +
  geom_boxplot(aes(x=category, y=allelic.div), alpha = 0, outlier.shape = " ") +
  geom_point( data=subset(df, allelic.div < .1 | allelic.div > 0.02),aes(x=category, y=allelic.div, colour=drug.resistance), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) + 
  geom_point( data=subset(df, allelic.div > .1 | allelic.div < 0.02),aes(x=category, y=allelic.div, colour=drug.resistance), size=5, alpha=.5) +
  coord_flip() + scale_fill_manual(values=c("#56B4E9","#E69F00")) + 
  guides(col = guide_legend(title = "Drug Resistance")) +
  theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$allelic.div), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by association with drug resistance") +
  xlab("")  + ylab("Alleles per nucleotidee") +
  geom_text( data=subset(df, allelic.div > .1 | allelic.div < 0.02), aes(x=category, y=allelic.div, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# with null model measure
ggplot(df) +
  geom_boxplot(aes(x=category, y=div.diff), alpha = 0, outlier.shape = " ") +
  geom_point( data=subset(df, div.diff < 0.09 | div.diff > -0.03),aes(x=category, y=div.diff), colour="coral2", size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) + 
  geom_point( data=subset(df, div.diff > 0.09 | div.diff < -0.03),aes(x=category, y=div.diff), colour="coral2", size=5, alpha=.5) +
  coord_flip() + theme_minimal() + guides(fill=FALSE) +
  theme(axis.text.y=element_text(size=14), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$div.diff), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category") +
  xlab("")  + ylab("Difference in alleles/nucleotide from expected based on mean") +
  geom_text( data=subset(df, div.diff > 0.09 | div.diff < -0.03), aes(x=category, y=div.diff, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by proportion of variable sites (nuc)
ggplot(df) +
  geom_boxplot(aes(x=category, y=p.vs.nuc), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=p.vs.nuc), colour="coral2", size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .3), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$p.vs.nuc), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Proportion of polymorphic sites to nonpolymorphic sites") +
  geom_text( data=subset(df, p.vs.nuc > 0.2 | p.vs.nuc < 0.01), aes(x=category, y=p.vs.nuc, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by ratio of unique amino acid sequences to unique nucleotide 
ggplot(df) +
  geom_boxplot(aes(x=category, y=r.count), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=r.count), colour="coral2", size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.1, .85), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$r.count), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Ratio of unique alleles in nucleotide to amino acid format") +
  geom_text( data=subset(df, r.count > .92 | r.count < .3), aes(x=category, y=r.count, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 

# by ratio of variable sites in nucleotide or amino acid form
ggplot(df) +
  geom_boxplot(aes(x=category, y=r.vs), alpha = 0, outlier.shape = " ") +
  geom_point( aes(x=category, y=r.vs, colour=log(avg.length)), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + scale_colour_gradientn(colours=rainbow(4)) + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$r.vs), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Ratio of variable sites in nucletodie to ami alleles in nucleotide to amino acid format") +
  geom_text( data=subset(df, r.vs > 2.8 | r.vs < 1), aes(x=category, y=r.vs, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 




ggplot(df, aes(p.vs.aa, allelic.div)) + geom_point() +
  geom_text( data=subset(df, p.vs.aa < 0.01), aes(x=p.vs.aa, y=allelic.div, label=name), size=4)

ggplot(df, aes(p.vs.nuc, div.diff)) + geom_point() +
  geom_text( data=subset(df, div.diff > 0.15), aes(x=p.vs.nuc, y=div.diff, label=name), size=4)

