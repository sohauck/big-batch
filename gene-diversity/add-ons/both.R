
filepath <- file.choose()

FIRSTheadtext <- scan(filepath, what = "complex", sep = "\n", n = 5)

# isolatesFIRST = as.numeric( strsplit(headtext[4], " ")[[1]][length(strsplit(headtext[4], " ")[[1]])] )

FIRST <- read.table(filepath, # open the file once it's in
           header=TRUE, sep="\t", # tab-separated with a header
           quote='"', skip = 5)


filepath <- file.choose()

SECONDheadtext <- scan(filepath, what = "complex", sep = "\n", n = 5)

# isolatesSECOND = as.numeric( strsplit(headtext[4], " ")[[1]][length(strsplit(headtext[4], " ")[[1]])] )

SECOND <- read.table(filepath, # open the file once it's in
                    header=TRUE, sep="\t", # tab-separated with a header
                    quote='"', skip = 5)


library(corrplot)

is.numeric(FIRST)

F.F <- subset(FIRST,  select = c(1:9) )
S.F <- subset(SECOND, select = c(1:9) )

B.F <- merge(F.F, S.F, by="Locus")
B.F <- subset(B.F, select = -c(1) )

M.F <- cor(B.F)
corrplot(M.F)

B <- merge(FIRST, SECOND, by="Locus")
B <- subset(B, select = -c(1) )

qplot(data = B.F, x = AvgLength.x, y = CountNuc.x)

M <- cor(B)

corrplot(B)
