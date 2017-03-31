library(ggplot2) # setting up packages
Args <- commandArgs(TRUE);      # retrieve args

setwd(Args[1]) # Going to directory with results
z <- ifelse(Args[2]=="TRUE",TRUE,FALSE) # turn to TRUE if want to plot by z-score instead
cat <- ifelse(Args[3]=="notdef",FALSE,TRUE) # TRUE if plotting by category
lbl <- ifelse(!is.na(Args[4]),as.numeric(Args[4]),10) # how many loci at both tail and head to label 


# Reading in table with output from Perl runs
df <- read.table("ResultsTable.txt", header = TRUE, sep = "\t")

# Adding locus categories if they are included
if ( cat == TRUE ) { 
  df.cat <- read.table(Args[3], header = TRUE, sep = "\t")
  names(df.cat)[1] <- "Locus"
  names(df.cat)[2] <- "Category"
  df <- merge( df, df.cat, by="Locus")
}

# Setting up new variables 
df$AllelicDiv <- c( df$CountNuc / df$AvgLength )
df$ADivNM <- c( (df$CountNuc - (mean(df$CountNuc / df$AvgLength) * df$AvgLength) ) / df$AvgLength )

df$VSitesNuc <- (df$AvgLength     - df$NonVarNuc) /  df$AvgLength
df$VSitesAA <- ((df$AvgLength / 3) - df$NonVarAA) / (df$AvgLength / 3)

df$RatioCount <- df$CountAA  /  df$CountNuc
df$RatioVS    <- df$VSitesAA /  df$VSitesNuc


# Writing fuller table back onto results directory
write.table (df, "CalculationsTable.txt", sep = "\t", row.names = FALSE)

# Moving to Graphs folder for saving graphs
setwd( paste (getwd(), "/Graphs/", sep = "") )

# =============
# Making graphs!
# =============

# Histogram of loci with >90% missing loci (actually all, but with dashed line)
ggplot( df, aes(x=Missing) ) + 
  geom_histogram( binwidth=(nrow(df)/30) ) +
  geom_vline( xintercept=(9/10)*nrow(df), 
              size=1, colour="red", linetype="dashed") + 
  theme_minimal() + 
  scale_x_continuous( limits=c( 0, nrow(df) ),
                      breaks=c( (1:4*(1/4) )*nrow(df) ),
                      labels=c("75%","50%","25%","00%")) +
  theme( axis.text.y=element_text(size=14), 
         axis.text.x=element_text(size=14),
         plot.title = element_text(face="bold"), 
         axis.title.x=element_text(vjust=-.5, size=14)) +
  ggtitle("Distribution of loci by percentage of isolates where they have a known sequence") +
    ylab("") + 
    xlab("Percentage of isolates in which locus is known")

ggsave("Dist-Missing.png", height = 9, width = 12, dpi = 100)

# Cutting off below 10% tagged 
df.all <- df # make back-up
cutoff <- (9/10) * nrow(df.all) # excluding anything tagged in less than 10% of isolate records
df <- df[df$Missing < cutoff,]
removed <- length(df.all[df.all$Missing > cutoff,1])
print ( paste ("Number of loci excluded due to being tagged in less than 10% of isolates:", removed) )


# Distribution of AllelicDiv and ADivNM

if ( z == TRUE ) {
  x1 <- scale(df$AllelicDiv, center = TRUE, scale = TRUE)
} else { x1 <- df$AllelicDiv }

ggplot( df, 
        aes( x = x1, 
             fill = "")) + 
  geom_histogram( alpha=.7, 
                  binwidth=((1/30)*(range(x1)[2]-range(x1)[1]))) + 
  ggtitle( "Density of alleles per nucleotide position" ) + 
  geom_vline( xintercept=mean(x1), 
              colour = "red", linetype="dashed") +
  theme_minimal( ) + 
  theme( plot.title = element_text(face="bold"), 
         axis.title.x=element_text(vjust=-.5, size=14)) +
  ylab("Density") + 
  xlab("Alleles per nucleotide letter") + 
  guides(fill=FALSE) 

ggsave("Dist-AllelicDiv.png", height = 9, width = 12, dpi = 100)


if ( z == TRUE ) {
  x1 <- scale(df$ADivNM, center = TRUE, scale = TRUE)
} else { x1 <- df$ADivNM }

ggplot( df, aes( x = x1, 
                 fill = "")) + 
  geom_histogram( alpha=.7, 
                  binwidth=((1/30)*(range(x1)[2]-range(x1)[1]))) + 
  geom_vline( xintercept = mean(x1), 
              colour = "red", linetype = "dashed") +
  ggtitle("Density of alleles per nucleotide position, compared to null model") +
  theme_minimal() + 
  theme( plot.title = element_text(face="bold"), 
         axis.title.x=element_text(vjust=-.5, size=14)) +
  ylab("Density") + 
  xlab("Difference in alleles per length in nucleotides from genome average") + 
  guides(fill=FALSE) 

ggsave("Dist-ADivNM.png", height = 9, width = 12, dpi = 100)


# The famous bubble plots!

# Point-AllelicDiv.png

if ( z == TRUE ) {
  y1 <- scale( df$AllelicDiv, center = TRUE, scale = TRUE)
} else { y1 <- df$AllelicDiv }

if ( cat == TRUE ) {
  x1 <- df$Category
} else { x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE) }

ggplot(df) +
  geom_point( data=df, 
              aes( x=x1, 
                   y=y1 ), 
              size=5, alpha=.5, 
              colour = "coral2",
              position = position_jitter(w=ifelse(cat,1/nlevels(df$Category),0) ) ) +
  geom_hline( yintercept = mean(y1), 
              size=1, colour="blue", linetype="dashed" ) +
  geom_text( aes( x=x1, 
                  y=y1, 
                  label=ifelse( y1 < head(sort(y1),lbl)[lbl] | 
                                y1 > tail(sort(y1),lbl)[1],
                                          as.character(Locus),'') ), 
             size=4, alpha=.8, vjust=-.5, angle = 30) +
  ggtitle("Genetic diversity of loci") + 
  xlab("") + 
  ylab("Alleles per nucleotide") +  
  coord_flip() + 
  theme_minimal() + 
  theme( axis.text.y = element_text(size=ifelse(cat,14,0)),
         axis.ticks  = element_line(size=ifelse(cat,0.5,0)),
         plot.title = element_text(face="bold"), 
         axis.title.x = element_text(vjust=-.5, size=14) )

ggsave("Point-AllelicDiv.png", height = 9, width = 12, dpi = 100)


# Point-ADivNM.png

if ( z == TRUE ) {
  y1 <- scale( df$ADivNM, center = TRUE, scale = TRUE)
} else { y1 <- df$ADivNM }

if ( cat == TRUE ) {
  x1 <- df$Category
} else { x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE) }

ggplot(df) +
  geom_point( data=df, 
              aes(x=x1, 
                  y=y1),
              size=5, alpha=.5,
              colour = "coral2",
              position = position_jitter(w=ifelse(cat,1/nlevels(df$Category),0) ) ) +
  geom_hline( yintercept = mean(y1),
              size=1, colour="blue", linetype="dashed") +
  geom_text( aes( x=x1,
                  y=y1, 
                  label=ifelse(   y1<head(sort(y1),lbl)[lbl] | 
                                  y1>tail(sort(y1),lbl)[1],
                                             as.character(Locus),'')), 
             size=4, alpha=.8, vjust=-.5, angle = 30) +
  coord_flip() + 
  ggtitle("Genetic diversity of loci, compared to null model") + 
  xlab("") + 
  ylab("Difference in alleles per nucleotide from genome average") +
  theme_minimal() + 
  theme( axis.text.y = element_text(size=ifelse(cat,14,0)),
         axis.ticks  = element_line(size=ifelse(cat,0.5,0)),
         plot.title = element_text(face="bold"), 
         axis.title.x = element_text(vjust=-.5, size=14) )

ggsave("Point-ADivNM.png", height = 9, width = 12, dpi = 100)


# Ratio of Count

if ( z == TRUE ) {
  y1 <- scale( df$RatioCount, center = TRUE, scale = TRUE)
} else { y1 <- df$RatioCount }

if ( cat == TRUE ) {
  x1 <- df$Category
} else { x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE) }

ggplot(df) +
  geom_point( data=df, 
              aes(x=x1, 
                  y=y1), 
              size=5, alpha=.5,
              colour = "coral2",
              position = position_jitter(w=ifelse(cat,1/nlevels(df$Category),0) ) ) +
  geom_hline( yintercept = mean(y1), 
              size=1, colour="blue", linetype="dashed") +
  geom_text( aes( x=x1, 
                  y=y1, 
                  label=ifelse( y1 < head(sort(y1),lbl)[lbl] | 
                                y1>= tail(sort(y1),lbl)[1],
                                             as.character(Locus),'')), 
             size=4, alpha=.8, vjust=-.5, angle = 30) +
  coord_flip() + 
  theme_minimal() + 
  theme( axis.text.y = element_text(size=ifelse(cat,14,0)), 
         plot.title = element_text(face="bold"), 
         axis.title.x=element_text(vjust=-.5, size=14), 
         axis.text.y = element_blank(),
         axis.ticks  = element_line(size=ifelse(cat,0.5,0))  ) +
  ggtitle("Distribution of purifying to diversifying selection") + 
  xlab("") + 
  ylab("Ratio of unique nucleotide to unique amino acid sequences per locus")

ggsave("Point-RatioCount.png", height = 9, width = 12, dpi = 100)


# Ratio of variable sites

if ( z == TRUE ) {
  y1 <- scale( df$RatioVS, center = TRUE, scale = TRUE)
} else { y1 <- df$RatioVS }

if ( cat == TRUE ) {
  x1 <- df$Category
} else { x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE) }
 
ggplot(df) +
  geom_point( data=df , 
              aes(x=x1, 
                  y=y1), 
              size=5, alpha=.5,
              colour = "coral2",
              position = position_jitter(w=ifelse(cat,1/nlevels(df$Category),0) ) ) +
  geom_hline( yintercept = mean(y1), 
             size=1, colour="blue", linetype="dashed") + 
  geom_text( aes( x=x1, 
                  y=y1, 
                  label=ifelse(   y1 < head(sort(y1),lbl)[lbl] | 
                                  y1 >=tail(sort(y1),lbl)[1],
                              as.character(Locus),'')), 
             size=4, alpha=.8, vjust=-.5, angle = 30) + 
  coord_flip() + 
  ggtitle("Distribution of purifying to diversifying selection") + 
  xlab("") + 
  ylab("Ratio of variable sites in nucleotide to amino acid format per locus") +
  theme_minimal() + 
  theme( axis.text.y = element_text(size=ifelse(cat,14,0)), 
         plot.title = element_text(face="bold"), 
         axis.title.x=element_text(vjust=-.5, size=14), 
         axis.text.y = element_blank(),
         axis.ticks  = element_line(size=ifelse(cat,0.5,0))  ) +

ggsave("Point-RatioVS.png", height = 9, width = 12, dpi = 100)
