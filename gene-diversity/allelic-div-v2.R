library(ggplot2)
library(grid)

# open table of isolates by loci filled with allele numbers
# no paralogous, 0 in for missing
GCtable <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Data sets/full.csv')
GCtable <- data.frame(GCtable[,-1], row.names=GCtable[,1])
GCtable[] <- lapply(GCtable, factor)

# creating table of allele counts
factorlevels <- vector(mode="numeric", length=0)
for (i in 1:ncol(GCtable) ) { 
  factorlevels <- c(factorlevels, nlevels(GCtable[,i]))
}
df <- as.data.frame( cbind(colnames(GCtable), as.numeric(factorlevels)))
names(df)[1] <- "locus"
names(df)[2] <- "allele.count"
df$allele.count <- as.numeric(levels(df$allele.count))[as.integer(df$allele.count)]

# bringing in more variables from lookup tables
dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/id-name-cat-conv.csv')
names(dfref)[1] <- "locus"
df <- merge(df,dfref,by="locus")
dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/locus-length.csv')
names(dfref)[2] <- "length"
df <- merge(df,dfref,by="locus")
dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Ref/sharing.txt')
df <- merge(df,dfref,by="locus")
dfref <- read.csv('/Volumes/sofia/Mycobacterium/MTBC allelic diversity/ECCMID/Data sets/full-missing.txt')
df <- merge(df,dfref,by="locus")

# setting up variables
df$allelic.div <- c( log10(  df$allele.count / df$length) )

# removing worst loci in terms of tagging (usually ones with seed problems)
cutoff <- 5000
df <- df[df$missing < cutoff,]

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

# for null model 
df$bef.div <- c( df$allele.count / df$length )
df$exp.div <- c( mean(df$bef.div) * df$length )
df$div.diff <- c( (df$allele.count - df$exp.div) / df$length )

ggplot(df, aes(x=div.diff)) + geom_histogram(alpha=.7, binwidth = 0.002) +
  geom_vline(xintercept=-mean(df$bef.div))

ggplot(df, aes(x=div.diff)) + stat_function(geom = "line", fun = dpois, arg = list(lambda=0), colour = "red")


ggplot(data.frame(X = seq(0, 30)), aes(x = X)) +
  stat_function(geom="line",fun=dpois,size=2, color="blue3", args =  list(lambda = 15), n = 31)



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
  geom_point( aes(x=category, y=div.diff, colour=drug.resistance), size=5, alpha=.5, position = position_jitter(w=0.15, h=0)) +
  coord_flip() + theme_minimal() + 
  theme(axis.text.y=element_text(size=14), legend.position=c(.9, .9), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  geom_hline(yintercept=mean(df$div.diff), size=1, colour="blue", linetype="dashed") +
  ggtitle("Genetic diversity of loci, split by locus functional category and coloured by gene length") +
  xlab("")  + ylab("Difference in alleles/nucleotide from excepted based on mean") +
  geom_text( data=subset(df, div.diff > 0.07 | div.diff < -0.03), aes(x=category, y=div.diff, label=name), size=4, alpha=.8, vjust=-1.5, position = position_jitter(width=0.3)) 






# Proportion of each category according to sharing
s.table <- aggregate(sharing2 ~ category, data = df, FUN = length)      

aggregate(x ~ Year + Month, data = df, FUN = length)

ggplot(df, aes(x=category, fill=sharing2)) + 
  geom_bar(position="fill", alpha=.7) + theme_minimal() +
  geom_hline(yintercept = .423, size=1) +
  scale_fill_manual(values=c("#56B4E9","#E69F00"), name = "Sharing status") +
  theme(axis.text.y = element_text(size=14), plot.title = element_text(face="bold"), axis.title.x=element_text(vjust=-.5, size=14)) +
  coord_flip(ylim = c(0, 1)) +
  ggtitle("Proportion of loci in each category with alleles specific to either side of deep branch") +
  xlab("")  + ylab("") +
  annotate("text", x = 1:11, y = 0.8, label = c(paste("p = ",res.fish[1,])), hjust = 0) +
  annotate("text", x = 1:11, y = 0.66, label = c("n.s.", "***","n.s.","n.s.","***","**","**","n.s.","*","n.s.","**"), hjust = 0)


# http://stats.stackexchange.com/questions/60764/cell-chi-square-test
# Partitioned chi-square test
M <- as.table(table(df$sharing2, df$category))

dimnames(M) <- list(sharing = c(levels(df$sharing2)),
                    category = c(levels(df$category)))

res.fish <- matrix(NA, nrow=dim(M)[1], ncol=dim(M)[2])
dimnames(res.fish) <- dimnames(M)

# Adapting for use of fisher's exact test instead, give odds ratio as result?
for ( i in 1:dim(M)[1] ) {
  for ( j in 1:dim(M)[2] ) {
    temp.table <- matrix(NA, 2, 2) # the collapsed 2x2 table
    
    temp.table[1,1] <- M[i,j]
    temp.table[1,2] <- sum(M[-i, j], na.rm=TRUE)
    temp.table[2,1] <- sum(M[i, -j], na.rm=TRUE)
    temp.table[2,2] <- sum(M[-i, -j], na.rm=TRUE)
    
    fish.t <- fisher.test(temp.table)
    
    res.fish[i, j] <- paste(round(fish.t$p.value, digits=3))  
  } 
} 




# Quantile-quantile plot, all and by cat.
qplot(sample = allelic.div, colour = factor(category), data = subset(df, category!="Unknown" & category!= "Ribosomal protein")) + 
  geom_abline(slope=.25, intercept=-3.05, alpha=.5)

indplot <- function(catnum, barcol) {
  p <- qplot(sample = allelic.div, data = subset(df, category == levels(df$category)[catnum]) ) + stat_qq(colour = barcol) +
    geom_abline(slope=.25, intercept=-3, alpha=.5) +
    xlab(levels(df$category)[catnum]) + ylab("") +
    theme(legend.position="none") + xlim(-3.5,3.5) + ylim(-5,-.5)
  return(p)
}

p1 <- indplot(1, "#F8766D")
p2 <- indplot(2, "#D39200")
p3 <- indplot(3, "#93AA00")
p4 <- indplot(4, "#00BA38")
p5 <- indplot(5, "#00C19F")
p6 <- indplot(6, "#00B9E3")
p7 <- indplot(7, "#619CFF")
p8 <- indplot(8, "#DB72FB")
p9 <- indplot(11, "#FF61C3")

multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, cols=3)



res.chi <- matrix(NA, nrow=dim(M)[1], ncol=dim(M)[2])
dimnames(res.chi) <- dimnames(M)

# Loop over all cells of the table, using chi-square, original
for ( i in 1:dim(M)[1] ) {
  for ( j in 1:dim(M)[2] ) {
    temp.table <- matrix(NA, 2, 2) # the collapsed 2x2 table
    
    temp.table[1,1] <- M[i,j]
    temp.table[1,2] <- sum(M[-i, j], na.rm=TRUE)
    temp.table[2,1] <- sum(M[i, -j], na.rm=TRUE)
    temp.table[2,2] <- sum(M[-i, -j], na.rm=TRUE)
    
    chi2 <- chisq.test(temp.table, correct=TRUE) # chi2-test with continuity correction
    
    # Automatically choose significance level (see SPSS documentation)
    sig.level <- ifelse(M[i,j] <= 300, 0.1,
                        ifelse(M[i,j] > 300 & M[i,j] <= 1000, 0.05,
                               ifelse(M[i,j] > 1000 & M[i,j] <= 4000, 0.025,
                                      ifelse(M[i,j] > 4000 & M[i,j] <= 20000,0.005,0.001)
                               ) ) )
    
    if ( chi2$p.value < sig.level ) {
      
      res.chi[i, j] <- paste(
        round(chi2$p.value, digits = 5), 
        chi2$observed[1],
        ifelse(chi2$observed[1] < chi2$expected[1], "<", ">"),
        round(chi2$expected[1], 2))
      
    } else {
      
      res.chi[i, j] <- paste( "n.s.", round(chi2$p.value, digits = 5))
      
    }
  } 
}
