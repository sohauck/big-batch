library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  # reading the table in to one whole data frame
  df.all <- reactive({
    inFile <- input$resultstable
    
    validate( need(input$resultstable != "", "Please upload a table on the left") )
    
    read.table(inFile$datapath, # open the file once it's in
                     header=TRUE, sep="\t", # tab-separated with a header
                     quote='"')
  })
  
  # subselecting the data frame depending on level of exclusion 
  df <- reactive({
    df.all <- df.all()
    
    cutoff <- ((100-input$percexc)/100)*nrow(df.all)
    df <- df.all[df.all$Missing < cutoff,]

    total.isolates <<- nrow(df.all)
    removed.isolates <<- nrow(df.all) - nrow(df)
    remain.isolates <<- total.isolates - removed.isolates
    
    return ( df )
  })

  output$exclusions <- renderText({
    
    total.isolates   <- nrow(df.all())
    removed.isolates <- nrow(df.all()) - nrow(df())
    remain.isolates <- (total.isolates - removed.isolates)
    
    paste(remain.isolates, 'loci;', removed.isolates, 'excluded out of', total.isolates) 
  })
  
  # need to change this for if there is a Category column in df
  # So that the "Use categories" checkbox only appears if a file is uploaded
  output$catinc <- reactive({
    return("Category" %in% colnames(df()))
  })
  outputOptions(output, 'catinc', suspendWhenHidden=FALSE)
  

  df.sel <- reactive ({
    # making visible copies of the relevant large data table
    df <- df()

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

    if ( input$categ == TRUE ) { df.sel$Category <- df$Category }
    
    return ( df.sel )
  })
  
  # if you want to check that the selected table looks like instead of the plot
  output$data.table <- renderDataTable({ df() })
  
  output$distplot <- renderPlot({
    
    # make a copy of the table that will be used
    df <- df.all()
    
    p <- ggplot( df, aes(x=Missing) ) +
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
      xlab("Distribution of number of isolates for which alocus has no known sequence")

    return ( p )
  })
  
  output$corrplot <- renderPlot({
    
    # make a copy of the table that will be used
    df <- df()
    
    qplot (Missing, AllelicDiv, data = df)
    
  })
  
  scatterPlot <- reactive({
    
    # make a copy of the table that will be used
    df.sel <- df.sel()

    p <- ggplot(df.sel) +
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
    
    if ( input$yVar == "AllelicDiv" ) 
      { p <- p + ylab("Alleles per nucleotide") }
    
    if ( input$yVar == "ADivNM" ) 
      { p <- p + ylab("Difference in alleles per length in nucleotides from genome average") }
    
    if ( input$yVar == "VSitesNuc" ) 
      { p <- p + ylab("Proportion of sites in nucleotide alignment which show any variation") }
    
    if ( input$yVar == "VSitesAA" ) 
    { p <- p + ylab("Proportion of sites in nucleotide alignment which show any variation") }
    
    if ( input$yVar == "RatioCount" )
      { p <- p + ylab("Ratio of unique nucleotide to unique amino acid sequences per locus") }
    
    if ( input$yVar == "RatioVS" ) 
      { p <- p + ylab("Ratio of variable sites in nucleotide to amino acid format per locus") }
    
    return ( p )
  })
  
  
  # what actually does the plot
  output$mainplot <- renderPlot({
    ggsave( "plot.pdf", scatterPlot() )
    return ( scatterPlot() )
  })
  
  output$info <- renderText({
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    
    paste0(
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "plot.pdf"
    },
    content = function(file) {
      file.copy("plot.pdf", file, overwrite=TRUE)
    }
  )
})

