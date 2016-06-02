df <- read.table(file.choose(), header=TRUE, sep="\t", quote='"', skip = 5)

df <- df[df$AllelicDiv < 0.2, ]

df.lm <- lm(AllelicDiv ~ Missing, data=df) 
df.lm.value <- summary(df.lm )$r.squared 

cutoff <- 8000

df.corr <- data.frame(x = numeric(800), y = numeric(800), stringsAsFactors = FALSE)
df.corr$x[1] <- cutoff
df.corr$y[1] <- df.lm.value


for ( i in 1 : 800 )
{
  cutoff <- cutoff - 10
  df.cut <- df[df$Missing < cutoff, ]
  df.corr$x[i] <- cutoff
  
  cut.lm <- lm(AllelicDiv ~ Missing, data=df.cut) 
  rs.value <- summary( cut.lm )$r.squared
  df.corr$y[i] <- rs.value
}


plot(  df.corr$y ~ df.corr$x )