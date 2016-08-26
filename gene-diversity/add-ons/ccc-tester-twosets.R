library(epiR)
library(ggplot2)

# import the big dataset
df1 <- read.table(file.choose(), header=TRUE, sep="\t", quote='"', skip = 5)
# remove Category if it is there
df1 <- subset(df1, select = -Category)

# take out excluded loci
df1 <- df1[df1$Missing < 7670-1143, ]


# make variables
df1$AllelicDiv <- c( df1$CountNuc / df1$AvgLength )
df1$ADivNM <- c( (df1$CountNuc - (mean(df1$CountNuc / df1$AvgLength) * df1$AvgLength) ) / df1$AvgLength )

df1$VSitesNuc <- (df1$AvgLength     - df1$NonVarNuc) /  df1$AvgLength
df1$VSitesAA <- ((df1$AvgLength / 3) - df1$NonVarAA) / (df1$AvgLength / 3)

df1$RatioCount <- df1$CountAA  /  df1$CountNuc
df1$RatioVS    <- df1$VSitesAA /  df1$VSitesNuc



df2 <- read.table(file.choose(), header=TRUE, sep="\t", quote='"', skip = 5)



df2$AllelicDiv <- c( df2$CountNuc / df2$AvgLength )
df2$ADivNM <- c( (df2$CountNuc - (mean(df2$CountNuc / df2$AvgLength) * df2$AvgLength) ) / df2$AvgLength )

df2$VSitesNuc <- (df2$AvgLength     - df2$NonVarNuc) /  df2$AvgLength
df2$VSitesAA <- ((df2$AvgLength / 3) - df2$NonVarAA) / (df2$AvgLength / 3)

df2$RatioCount <- df2$CountAA  /  df2$CountNuc
df2$RatioVS    <- df2$VSitesAA /  df2$VSitesNuc

# subset based on loci that survived major cut

df3 <- merge(df1, df2, by = "Locus")

res <- epi.ccc(df3$AllelicDiv.x, df3$AllelicDiv.y, ci = "z-transform", conf.level = 0.95)

qplot(df3$AllelicDiv.x, df3$AllelicDiv.y)

for ( i in 1 : length(names(df2)))
{
  epi.ccc(df1[i], df2[i], ci = "z-transform", conf.level = 0.95)
}

