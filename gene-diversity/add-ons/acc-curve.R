library(ggplot2)

path = "/Users/sofia/GitHub/big-batch/gene-diversity/add-ons/r-s-c"

file.names <- dir(path, pattern =".txt", full.names = TRUE)

vectorx <- vector()
vectorx[1] = 0

for(i in 1:length(file.names)) {
  df <- read.table(  file = paste ("/Users/sofia/GitHub/big-batch/gene-diversity/add-ons/r-s-c/Res-CountNuc-accumulate-", i, ".txt.txt" , sep="" ), # open the file once it's in
                     header=TRUE, sep="\t", # tab-separated with a header
                     quote='"', skip = 3, row.names = NULL)
  
  vectorx[i+1] <- sum( df[,4] )
}


qplot(seq_along(vectorx), vectorx)
