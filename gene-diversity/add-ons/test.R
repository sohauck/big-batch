path = "/Volumes/sofia/Manuscripts/CompBio/Mtub-Data/acc20-CountNucs-5"

file.names <- dir(path, pattern =".txt", full.names = TRUE)

vectorx <- vector()
vectorx[1] = 0

for(i in 1:length(file.names)) {
  #file <- read.table(file.names[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
  df <- read.table(file.names[i], # open the file once it's in
                     header=TRUE, sep="\t", # tab-separated with a header
                     quote='"', skip = 3)
  
  #out.file <- rbind(out.file, file)
  vectorx[i+1] <- vectorx[i] + sum( df$CountNuc )
}

write.table(out.file, file = "cand_Brazil.txt",sep=";", 
            row.names = FALSE, qmethod = "double",fileEncoding="windows-1252")


df2 <- read.table(file.choose(), header = TRUE)
                 
                 
                 # open the file once it's in
                 header=TRUE, sep="\t", # tab-separated with a header
                 quote='"', skip = 3)
