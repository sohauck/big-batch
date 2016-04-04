library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  df <- reactive({
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
    
    return(df)
  })

  df.cat <- reactive({
    inFile2 <- input$categorytable

    if (is.null(inFile2)) # until have a file,
    { return(NULL) }
    else {
      df.cat <- read.table(inFile2$datapath,
                           header = TRUE,
                           sep = "\t")

      names(df.cat)[1] <- "Locus"
      names(df.cat)[2] <- "Category"
      
      return(df.cat)
    }
  })

  output$catfileUp <- reactive({
    return(!is.null(input$categorytable))
  })

  outputOptions(output, 'catfileUp', suspendWhenHidden=FALSE)

  df.sel <- reactive ({
    
    # making visible copies of the relevant tables
    df     <- df()
    
    #
    if (!is.null(input$categorytable)) { # until have a file,
      df.cat <- df.cat()
      df <- merge( df, df.cat, by = "Locus")
    }
    
    if ( input$categ == TRUE ) {
      x1 <- df$Category }
    else { x1 <- sample(seq(from=0, to=1, by=.01), size = 25, replace = TRUE) }
    
      if ( input$yVar == "AllelicDiv" ) {
        if ( input$zscore == FALSE ) { y1 <- df$AllelicDiv }
        else { y1 <- scale( df$AllelicDiv, center = TRUE, scale = TRUE) }
        }

      if ( input$yVar == "ADivNM" ) {
        if ( input$zscore == FALSE ) { y1 <- df$ADivNM }
        else { y1 <- scale( df$ADivNM, center = TRUE, scale = TRUE) }
      }

      if ( input$yVar == "VSitesNuc" ) {
        if ( input$zscore == FALSE ) { y1 <- df$VSitesNuc }
        else { y1 <- scale( df$VSitesNuc, center = TRUE, scale = TRUE) }
      }

      if ( input$yVar == "VSitesAA" ) {
        if ( input$zscore == FALSE ) { y1 <- df$VSitesAA }
        else { y1 <- scale( df$VSitesAA, center = TRUE, scale = TRUE) }
      }

      if ( input$yVar == "RatioCount" ) {
        if ( input$zscore == FALSE ) { y1 <- df$RatioCount }
        else { y1 <- scale( df$RatioCount, center = TRUE, scale = TRUE) }
      }

      if ( input$yVar == "RatioVS" ) {
        if ( input$zscore == FALSE ) { y1 <- df$RatioVS }
        else { y1 <- scale( df$RatioVS, center = TRUE, scale = TRUE) }
      }
    
    Locus <- df$Locus
    

    df.sel <- data.frame( Locus, x1, y1 )
    colnames(df.sel) <- c( "Locus", "xsel", "ysel" )

    
    return ( df.sel )

  })
  
#   output$my_output_data <- renderTable({ df.sel() })  

  output$mainplot <- renderPlot({

    df.sel <- df.sel()

    # some way to plot distributions?

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

    # ifelse(input$categ,factor(xsel),

    ggplot(df.sel) +
      geom_point( data=df.sel,
                  aes( x=xsel,
                       y=ysel
                       #colour = ifelse(use.cat,factor(xsel),"coral2")
                       ),
                  size=5, alpha=.5,
                  position = position_jitter(w=.5) ) +
      geom_hline( yintercept = mean(df.sel$ysel),
                  size=1, colour="blue", linetype="dashed" ) +
      ggtitle("Genetic diversity of loci") +
      xlab("") +
      ylab("Alleles per nucleotide") +
      coord_flip() +
      theme_minimal() +
      theme( axis.text.y = element_text(size=ifelse(input$categ,14,0)),
             axis.ticks  = element_line(size=ifelse(input$categ,0.5,0)),
             plot.title = element_text(face="bold"),
             axis.title.x = element_text(vjust=-.5, size=14) )
    # geom_text( aes( x=x1,
    #                 y=y1,
    #                 label=ifelse( y1 < head(sort(y1),lbl)[lbl] |
    #                                 y1 > tail(sort(y1),lbl)[1],
    #                               as.character(Locus),'') ),
    #            size=4, alpha=.8, vjust=-.5, angle = 30) +


  })
  
   
})

