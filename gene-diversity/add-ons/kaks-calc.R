library(seqinr)


alignment <- read.alignment(file=file.choose(),format="FASTA")
res <- kaks(alignment)
ma <- mean(res$ka-res$ks)
