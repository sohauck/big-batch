library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  output$mainplot <- renderPlot({
   
    inFile <- input$resultstable
    
    if (is.null(inFile)) # until have a file, 
      return(NULL) # just keep it all blank
    
    df <- read.table(inFile$datapath, # open the file once it's in
                         header=TRUE, sep="\t", # tab-separated with a header
                         quote='"')
    
    # Setting up new variables 
    df$AllelicDiv <- c( df$CountNuc / df$AvgLength )
    df$ADivNM <- c( (df$CountNuc - (mean(df$CountNuc / df$AvgLength) * df$AvgLength) ) / df$AvgLength )
    
    df$VSitesNuc <- (df$AvgLength     - df$NonVarNuc) /  df$AvgLength
    df$VSitesAA <- ((df$AvgLength / 3) - df$NonVarAA) / (df$AvgLength / 3)
    
    df$RatioCount <- df$CountAA  /  df$CountNuc
    df$RatioVS    <- df$VSitesAA /  df$VSitesNuc
    

    
    if ( input$yVar == "AllelicDiv" ) {
      if ( input$zscore == FALSE ) { y1 <- df$AllelicDiv }
      else { y1 <- scale( df$AllelicDiv, center = TRUE, scale = TRUE) }
      }
      
    if ( input$yVar == "ADivNM" ) {
      if ( input$zscore == FALSE ) { y1 <- df$ADivNM }
      else { y1 <- scale( df$ADivNM, center = TRUE, scale = TRUE) }
    }
    
    qplot(y1)
    
    # ggplot( df, aes(x=Missing) ) + 
    #   geom_histogram( binwidth=(nrow(df)/30) ) +
    #   geom_vline( xintercept=(9/10)*nrow(df), 
    #               size=1, colour="red", linetype="dashed") + 
    #   theme_minimal() + 
    #   scale_x_continuous( limits=c( 0, nrow(df) ),
    #                       breaks=c( (1:4*(1/4) )*nrow(df) ),
    #                       labels=c("75%","50%","25%","00%")) +
    #   theme( axis.text.y=element_text(size=14), 
    #          axis.text.x=element_text(size=14),
    #          plot.title = element_text(face="bold"), 
    #          axis.title.x=element_text(vjust=-.5, size=14)) +
    #   ggtitle("Distribution of loci by percentage of isolates where they have a known sequence") +
    #   ylab("") + 
    #   xlab("Percentage of isolates in which locus is known")
    
    
    # if ( z == TRUE ) {
    #   y1 <- scale( df$AllelicDiv, center = TRUE, scale = TRUE)
    # } else { y1 <- df$AllelicDiv }
    # 
    # if ( cat == TRUE ) {
    #   x1 <- df$Category
    # } else { 
    
    # x1 <- sample(seq(from = 0, to = 1, by = .01), size = nrow(df), replace = TRUE) 
    # 
    # 
    # ggplot(df) +
    #   geom_point( data=df, 
    #               aes( x=x1, 
    #                    y=AllelicDiv ), 
    #               size=5, alpha=.5, 
    #               colour = "coral2" ) +
    #   geom_hline( yintercept = mean(AllelicDiv), 
    #               size=1, colour="blue", linetype="dashed" ) +
    #   # geom_text( aes( x=x1, 
    #   #                 y=y1, 
    #   #                 label=ifelse( y1 < head(sort(y1),lbl)[lbl] | 
    #   #                                 y1 > tail(sort(y1),lbl)[1],
    #   #                               as.character(Locus),'') ), 
    #   #            size=4, alpha=.8, vjust=-.5, angle = 30) +
    #   ggtitle("Genetic diversity of loci") + 
    #   xlab("") + 
    #   ylab("Alleles per nucleotide") +  
    #   coord_flip() + 
    #   theme_minimal()
    #   # theme( axis.text.y = element_text(size=ifelse(cat,14,0)),
    #   #        axis.ticks  = element_line(size=ifelse(cat,0.5,0)),
    #   #        plot.title = element_text(face="bold"), 
    #   #        axis.title.x = element_text(vjust=-.5, size=14) )
    # 
    
  })
   
})

