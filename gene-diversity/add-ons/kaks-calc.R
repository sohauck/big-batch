library(seqinr)
library(corrplot)


# read full size
df1 <- read.table(file.choose(), header=TRUE, sep="\t", quote='"', skip = 5)

# take out excluded loci
df1 <- df1[df1$Missing < 7670-1143, ] #specific to RaMi dataset!!

# make variables
df1$AllelicDiv <- c( df1$CountNuc / df1$AvgLength )
df1$ADivNM <- c( (df1$CountNuc - (mean(df1$CountNuc / df1$AvgLength) * df1$AvgLength) ) / df1$AvgLength )

df1$VSitesNuc <- (df1$AvgLength     - df1$NonVarNuc) /  df1$AvgLength
df1$VSitesAA <- ((df1$AvgLength / 3) - df1$NonVarAA) / (df1$AvgLength / 3)

df1$RatioCount <- df1$CountAA  /  df1$CountNuc
df1$RatioVS    <- df1$VSitesAA /  df1$VSitesNuc


# find the folder with reverse aligned FASTA files
FASTAdir <- dirname(file.choose())

# get all FASTA files from folder
FASTAs <- list.files( path = FASTAdir, pattern = ".FAS", full.names = TRUE, ignore.case = TRUE)

# empty data frame for results
res.df <- data.frame()

for ( i in 1 : length (FASTAs) ) #for each file
{
  # open aligment
  temp.align <- read.alignment( file = FASTAs[i], format = "FASTA")

  # kaks it
  temp.kaks.list <- kaks(temp.align)
  
  # put locus name on the dataframe
  locusname <- strsplit( basename(FASTAs[i]) , '[.]')[[1]]
  res.df[i,1] <- locusname[1]
  
  # converting data frames to matrices for speed
  ma.ka <- as.matrix(temp.kaks.list$ka)
  ma.ks <- as.matrix(temp.kaks.list$ks)

  # correct negative values to 0 in matrix (since neg. means not calc, see kaks manual)
  ma.ka[ma.ka < 0] <- 0
  ma.ks[ma.ks < 0] <- 0
  
  # pre-calculations for speed
  sum.ma.ka <- sum(ma.ka)
  sum.ma.ks <- sum(ma.ks)
  
  # putting results into res.df frame
  res.df[i,2] <- sum.ma.ka + sum.ma.ks
  res.df[i,3] <- mean(ma.ka + ma.ks)

  res.df[i,4] <- sum.ma.ka - sum.ma.ks
  res.df[i,5] <- mean(ma.ka - ma.ks)

  res.df[i,6] <- sum.ma.ka / sum.ma.ks
  res.df[i,7] <- mean(ma.ka / ma.ks)

  res.df[i,8] <- sum.ma.ka 
  res.df[i,9] <- sum.ma.ks

  res.df[i,10] <- mean(ma.ka)
  res.df[i,11] <- mean(ma.ks)
  
  align.length <- nchar(temp.align$seq[[1]])
  
  res.df[i,12] <- res.df[i,2] / align.length
  res.df[i,13] <- res.df[i,4] / align.length
  res.df[i,14] <- res.df[i,6] / align.length
}

colnames(res.df) <- c("Locus",
                      "SumAdded","IndAdded","SumMinus","IndMinus","SumDivid","IndDivid",
                      "SumKS","SumKA","MeanKS","MeanKA","SumAddedbyL","SumMinusbyL","SumDividbyL")



# merge data sets
merged.df <- merge(df1, res.df, by = "Locus")


# removing infinite values
merged.df <- merged.df[-c(3772, 635), ] #MYCO010959, MYCO000657, KA = 0


 
# remove non-numeric for correlation testing
merged.df <- subset(merged.df, select = - c(Locus,IndDivid))
# removing IndDivid since all NaN

# making correlation plots
M <- cor(merged.df)
merged.df$LSumDividbyL <- log10(merged.df$SumDividbyL)
corrplot(M)

df.red1 <- subset(merged.df, select = c(RatioCount, SumDivid))
df.red1$LogRC <- log(df.red1$RatioCount)
df.red1$LogSD <- log(df.red1$SumDivid)

M.1 <- cor( df.red1 )
corrplot(M.1)


df.red2 <- subset(merged.df, select = c(AllelicDiv, SumAddedbyL))
df.red2$LogAD <- log(df.red2$AllelicDiv)
df.red2$LogSA <- log(df.red2$SumAddedbyL)

M.2 <- cor( df.red2 )
corrplot(M.2)


# AllelicDiv as measure of SumAddedbyL

# simple linear regression
qplot   (log(merged.df$AllelicDiv),
         log(merged.df$SumAddedbyL),
         alpha=merged.df$VSitesNuc)
cor.test(log(merged.df$AllelicDiv),log(merged.df$SumAddedbyL))
#Decent, R = 0.7763702 ,   0.7635075 0.7886171 at 95%

# simple linear regression
df.AD <- subset(merged.df, select = c(SumAddedbyL, AllelicDiv, VSitesNuc))

fit.AD <- lm( log(SumAddedbyL) ~ log(AllelicDiv) , data = df.AD ) # 0.6028
summary(fit.AD)

fit.AD2 <- lm( log(SumAddedbyL) ~ log(AllelicDiv) + (VSitesNuc), data = df.AD ) # 0.7755
summary(fit.AD2)


ggplot (data = merged.df, aes(x = log(merged.df$AllelicDiv), y = log(merged.df$SumAddedbyL) )) +
         geom_point() + geom_smooth()



# RatioCount as measure of SumDivid

# Plain
qplot   (merged.df$RatioCount,merged.df$SumDivid)

# Straighter
qplot   (merged.df$RatioCount,log(merged.df$SumDivid))
qplot   (merged.df$RatioCount,log(merged.df$SumDivid),colour=merged.df$VSitesNuc)
ggplot (data = merged.df, aes(x = merged.df$RatioCount, y = log(merged.df$SumDivid))
        ,colour=merged.df$VSitesNuc) + geom_point() + geom_smooth()

# Plain 
cor.test( merged.df$RatioCount , merged.df$SumDivid)  # 0.5552726

# Log transformations
cor.test(      merged.df$RatioCount, log10(merged.df$SumDivid)) # 0.8040767


df.RC <- subset(merged.df, select = c(SumDivid, RatioCount, AvgLength))

fit.RC <- lm( log(SumDivid) ~ RatioCount , data = df.RC ) # 0.6465,	Adjusted R-squared:  0.6464 
summary(fit.RC)
anova(fit.RC)

qplot(seq_along(fit.RC$residuals), fit.RC$residuals)

ggplot (data = merged.df, aes(x = merged.df$RatioCount, y = log(merged.df$SumDivid))) + 
  geom_point( colour = df.RC$AvgLength) + 
  stat_smooth( method = "lm" )#, formula = log(merged.df$SumDivid) ~ merged.df$RatioCount)
  geom_abline( intercept = -2.94958, slope = 4.08015 )

  
  qplot   (merged.df$RatioCount,
           log(merged.df$SumDivid),
           colour=merged.df$CountAA)
  
  
  

write.table(res.df, file = paste(dirname,"res-df.txt",sep = "") )
