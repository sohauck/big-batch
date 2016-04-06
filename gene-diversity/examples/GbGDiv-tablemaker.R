library(ggplot2) # setting up packages
Args <- commandArgs(TRUE);      # retrieve args

setwd(Args[1]) # Going to directory with results
cat <- ifelse(Args[2]=="notdef",FALSE,TRUE) # TRUE if plotting by category

# Reading in table with output from Perl runs
df <- read.table("ResultsTable.txt", header = TRUE, sep = "\t")


# Setting up new variables 
df$AllelicDiv <- c( df$CountNuc / df$AvgLength )
df$ADivNM <- c( (df$CountNuc - (mean(df$CountNuc / df$AvgLength) * df$AvgLength) ) / df$AvgLength )

df$VSitesNuc <- (df$AvgLength     - df$NonVarNuc) /  df$AvgLength
df$VSitesAA <- ((df$AvgLength / 3) - df$NonVarAA) / (df$AvgLength / 3)

df$RatioCount <- df$CountAA  /  df$CountNuc
df$RatioVS    <- df$VSitesAA /  df$VSitesNuc


# Adding locus categories if they are included
if ( cat == TRUE ) { 
  df.cat <- read.table(Args[2], header = TRUE, sep = "\t")
  names(df.cat)[1] <- "Locus"
  names(df.cat)[2] <- "Category"
  df <- merge( df, df.cat, by="Locus")
}


# Writing fuller table back onto results directory
write.table (df, "ResultsTable.txt", sep = "\t", row.names = FALSE)
