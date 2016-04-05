library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  df.all <- reactive({
    inFile <- input$resultstable
    
    validate(
      need(input$resultstable != "", "Please upload a table on the left")
    )
    
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
    
    # df.all <- df # make back-up
    # cutoff <- (9/10) * nrow(df.all) # excluding anything tagged in less than 10% of isolate records
    # df <- df[df$Missing < cutoff,]
    # removed <- length(df.all[df.all$Missing > cutoff,1])
    # print ( paste ("Number of loci excluded due to being tagged in less than 10% of isolates:", removed) )
    
    
    
    return( df )
  })
  
  df <- reactive({
    df.all <- df.all()
    
    cutoff <- ((100-input$percexc)/100)*nrow(df.all)
    df <- df.all[df.all$Missing < cutoff,]

    total.isolates <<- nrow(df.all)
    removed.isolates <<- nrow(df.all) - nrow(df)
    remain.isolates <<- total.isolates - removed.isolates
    
    return ( df )
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

  # So that the "Use categories" checkbox only appears if a file is uploaded
  output$catfileUp <- reactive({
    return(!is.null(input$categorytable))
  })
  outputOptions(output, 'catfileUp', suspendWhenHidden=FALSE)

  df.sel <- reactive ({
    # making visible copies of the relevant large data table
    df <- df()
    # and likewise for the category table if it exists, and then also merge it into the main
    if (!is.null(input$categorytable)) { # until have a file,
      df.cat <- df.cat()
      df <- merge( df, df.cat, by = "Locus")
    }
    
    # don't complain yet if nothing to compare to... 
    # if  (is.null(input$categorytable))
    # { return (NULL) }
    
    # make the x-variable into either the category or a random number between 0 and 1 for visibility
    if ( input$categ == TRUE ) {
      x1 <- df$Category }
    else { x1 <- sample(seq(from=0, to=1, by=.01), size = nrow(df), replace = TRUE) }
    
    # make the y-variable whatever the selected variable is, with z-score scaling if necessary 
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
    
    # copy the Locus column to a vector for merging into the selected data frame
    Locus <- df$Locus
    
    # make the selected data frame
    df.sel <- data.frame( Locus, x1, y1 )
    colnames(df.sel) <- c( "Locus", "xsel", "ysel" )

    return ( df.sel )
  })
  
  # if you want to check that the selected table looks like instead of the plot
  output$data.table <- renderTable({ df() })
  
  output$distplot <- renderPlot({
    
    # make a copy of the table that will be used
    df <- df.all()
    
    ggplot( df, aes(x=Missing) ) +
      geom_histogram( binwidth=(nrow(df)/30) ) +
      geom_vline( xintercept=((100-input$percexc)/100)*nrow(df),
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
  
  })
  
  output$corrplot <- renderPlot({
    
    # make a copy of the table that will be used
    df <- df()
    
    qplot (Missing, paste(input$yVar), data = df)
    
  })
  
  # what actually does the plot
  output$mainplot <- renderPlot({
  
    # make a copy of the table that will be used
    df.sel <- df.sel()

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

