first.file <- file.choose()

# directory where all the accumulated count files are
directory <- dirname(first.file)

# get all files from folder
allfiles <- list.files( path = directory, pattern = ".txt", full.names = TRUE, ignore.case = TRUE)


# start data frame with that first file
total.df <- read.table(first.file, sep="\t", header = TRUE, quote='"', skip = 2)
total.df <- subset(total.df, select = c("Locus", "CountNuc"))



# limit to wanted loci
df1 <- read.table(file.choose(), sep="\t", header = TRUE, quote='"', skip = 5)
df1 <- df1[df1$Missing < 7670-1143, ]

total.df <- total.df[total.df$Locus %in% df1$Locus, ]
df1.lim <- subset(df1, select = c("Locus", "AvgLength"))
total.df <- merge(total.df, df1.lim, by = "Locus")


total <- length(allfiles) 

for ( i in 2 : total ) #for each file
{
  this.file <- paste( substr(first.file, 1, nchar(first.file)-5), i, ".txt", sep = "")
  temp.df <- read.table(this.file, sep="\t", header = TRUE, quote='"', skip = 2)
  temp.df <- subset(temp.df, select = c("Locus", "CountNuc"))
  colnames(temp.df) <- c("Locus", i)
  total.df <- merge(temp.df, total.df, by = "Locus" )
  total.df[2] <- total.df[2] / total.df$AvgLength
  total.df[3] <- total.df[3] / total.df$AvgLength
}

total.df.reserve <- total.df

total.ma <- as.matrix(total.df)
div.ma <- total.ma %o% 1/df1$AvgLength

ccc.df <- data.frame()


for ( i in 2 : total ) #for each file
{
  zscale <- scale( total.df[,i], center = TRUE, scale = TRUE)
  temp.res <- epi.ccc(zscale, total.df[,total], ci = "z-transform", conf.level = 0.95)
  ccc.df[i,1] <- i
  temp.res.df <- temp.res$rho.c
  ccc.df[i,2] <- temp.res.df$est
  ccc.df[i,3] <- temp.res$s.shift
  ccc.df[i,4] <- temp.res$l.shift
  ccc.df[i,5] <- temp.res$C.b
  
  temp.res2 <- cor.test(total.df[,i], total.df[,total], conf.level = 0.95) #method = pearson, kendall, spearman
  temp.res.df <- temp.res2$estimate
  ccc.df[i,6] <- temp.res.df
  
  colnames(ccc.df) <- c("counter","Rho","Scale.s","Loc.s","C.b","Rest")
}

qplot( seq(from = 1, to = 50, by = 1), ccc.df$Rho )
qplot( seq(from = 1, to = 50, by = 1), ccc.df$Scale.s )
qplot( seq(from = 1, to = 50, by = 1), ccc.df$Loc.s )
qplot( seq(from = 1, to = 50, by = 1), ccc.df$C.b )
qplot( seq(from = 1, to = 50, by = 1), ccc.df$Rest )



qplot(seq_along(sort(ccc.df$Rho, decreasing = FALSE)), sort(ccc.df$Rho, decreasing = FALSE) )
qplot(seq_along(sort(ccc.df$Scale.s, decreasing = FALSE)), sort(ccc.df$Scale.s, decreasing = FALSE) )
qplot( seq(from = 1, to = 50, by = 1), ccc.df2$Scale.s )

qplot(seq_along(sort(ccc.df$Loc.s, decreasing = FALSE)), sort(ccc.df$Loc.s, decreasing = FALSE) )

qplot(ccc.df$V1, ccc.df$V5)

qplot(total.df[,4], total.df[,total])

