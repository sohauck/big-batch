library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  # reading the table in to one whole data frame
  df.all <- reactive({
    inFile <- input$resultstable
    
    validate( need(input$resultstable != "", "Please upload a table on the left") )
    
    headtext <- scan(inFile$datapath, what = "complex", sep = "\n", n = 5)
    
    isocount$total = as.numeric( strsplit(headtext[4], " ")[[1]][length(strsplit(headtext[4], " ")[[1]])] )
    
    read.table(inFile$datapath, # open the file once it's in
                     header=TRUE, sep="\t", # tab-separated with a header
                     quote='"', skip = 5)
  })
  
  # subselecting the data frame depending on level of exclusion 
  df <- reactive({
    df.all <- df.all()
    
    locount$total <- nrow(df.all)
    
    if ( input$percexc == 0 )
    { return ( df.all ) }
    else {
      # 
      isocount$cutoff <- floor((1-(input$percexc/100))*isocount$total)
      
      # excluse from data frame all those which are worse off that the cutoff value
      df <- df.all[df.all$Missing < isocount$cutoff,]
      
      validate( need(nrow(df) != "0", "You removed all loci. Please adjust the slider on the left.") )
      
      return ( df )
    }
  })

  # text output that describes total number of loci in data table and number excluded / remaining
  
  locount <- reactiveValues( total = 0, removed = 0, remain = 0)
  isocount <- reactiveValues( total = 0, cutoff = 0)
  
  output$exclusions <- renderText({
      locount$removed <- locount$total - nrow(df())
      locount$remain  <- locount$total - locount$removed
      
      line1 <- paste('Isolates reviewed:', isocount$total)
      line2 <- paste('Loci included:', locount$remain)
      line3 <- paste(locount$removed ,'out of possible', locount$total, 'loci were removed. Cutoff of', 
      input$percexc, '% maximum isolates with missing allele designation filtered out loci in less than', isocount$total - isocount$cutoff ,'isolates.')
      
      paste (line1, line2, line3, sep="\n")
  })
  
  # checks if Category column exists in table and creates "Use categories" checkbox if is so
  output$catinc <- reactive({
    return("Category" %in% colnames(df()))
  })
  outputOptions(output, 'catinc', suspendWhenHidden=FALSE)

  # makes little data table that will be used for the graph, taking into account categ, yVar and zscore
  df.sel <- reactive ({
    # making visible copies of the relevant large data table
    df <- df()

    # make the x-variable into either the category or a random number between 0 and 1 for visibility
    if ( input$categ == TRUE ) {
      x1 <- df$Category }
    else { x1 <- sample(seq(from=0, to=1, by=.001), size = nrow(df), replace = TRUE) }
    
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
    
    # make the selected data frame
    df.sel <- data.frame( df$Locus, df$Missing, df$Paralogous, df$AvgLength, x1, y1) #
    colnames(df.sel) <- c( "Locus", "Missing", "Paralogous", "AvgLength", "xsel", "ysel") # 

    if ( input$categ == TRUE ) { df.sel$Category <- df$Category }
    
    return ( df.sel )
  })
  
  # output of the data table, with exclusions already removed, and labelled points on column
  output$data.table <- renderDataTable({ 
    df <- df()
    df$Label<- "No"
    df$Label[df$Locus %in% labeltag$list] <- "Yes"
    return ( df )
    })
  
  # plots distribution of Missing values across entire dataset 
  output$distplot <- renderPlot({
    
    # make a copy of the table that will be used
    df <- df.all()
    
    p <- ggplot( df, aes(x=Missing) ) +
      geom_histogram( binwidth=(isocount$total/40) ) + # so that there are always 30 bins 
      geom_vline( xintercept=ifelse(input$percexc,isocount$cutoff,isocount$total),
                  size=1, colour="red", linetype="dashed") +
      coord_cartesian( xlim = c( 0, isocount$total ) ) +
      theme_minimal() + 
      ggtitle("Distribution of loci by number of isolates for which there was has no known sequence")
    
    return ( p )
  })
  
  # plots Missing vs Allelic Div 
  output$corrplot <- renderPlot({
    df <- df()
    
    ggplot( df, aes (x = Missing, y = AllelicDiv) )+
      geom_point() + 
      geom_smooth( span = .01*locount$total ) + # so that span scales with n 
      coord_cartesian( xlim = c( 0, isocount$total ) ) +
      theme_minimal() +
      ggtitle("Regression of 'missing' designation and allelic diversity value per locus")
     
  })
  
  brange <- reactiveValues( xmin = NULL, ymin = NULL, xmax = NULL, ymax = NULL )
  labeltag <- reactiveValues ( list = NULL, use = FALSE )
  
  observeEvent(input$resetlabel, {
    labeltag$list <- NULL
    labeltag$use <- FALSE 
  })
  
  # if a double click happens, add the brushed points' Locus ID to labeltag$list
  observeEvent( { input$mp_dblclick }, {
    
    brush <- input$mp_brush
    
    if (!is.null(brush)) {
      brange$xmin <- brush$xmin
      brange$ymin <- brush$ymin
      brange$xmax <- brush$xmax
      brange$ymax <- brush$ymax
    } else {
      brange$xmin <- NULL
      brange$ymin <- NULL
      brange$xmax <- NULL
      brange$ymax <- NULL
    }
    
    df.sel <- df.sel()
    
    # revise the labeltag list into non-duplicated of appended
    if (!input$categ) {
      templist <- append ( labeltag$list, 
                           as.character( subset (df.sel,
                                                 ysel  >  brange$xmin &
                                                 ysel  <  brange$xmax &
                                                 xsel  >  brange$ymin &
                                                 xsel  <  brange$ymax)[,"Locus"])) }
    
    else {
      templist <- append ( labeltag$list, 
                           as.character( subset ( df.sel, 
                                Category %in% levels(df.sel$Category)[round(brange$ymin):round(brange$ymax)] &
                                                  ysel  >  brange$xmin &
                                                  ysel  <  brange$xmax )[,"Locus"])) }
    
    labeltag$list <- templist[!duplicated(templist)]
    templist <- NULL
    
    # turn on the added label layer only if there are points to label, avoids error messages
    if (!is.null(labeltag$list))
    { labeltag$use <- TRUE }
  })
  
  # the main plot of the whole thing 
  scatterPlot <- reactive({
    
    # make a copy of the table that will be used
    df.sel <- df.sel()

    # add the colour variable 
    if ( input$colourVar == "None (labels only)" ) { 
      df.sel$csel <- ifelse(df.sel$Locus %in% labeltag$list, "Selected", "Non-selected")
      colourtitle <- ""
    } #
    
    if ( input$colourVar == "Missing" ) { 
      df.sel$csel <- df.sel$Missing / isocount$total 
      colourtitle <- "Percentage of isolates where locus is not tagged"
    }
    
    if ( input$colourVar == "Paralogous" )
    { 
      df.sel$csel <- df.sel$Paralogous / isocount$total
      colourtitle <- "Percentage of isolates where locus is paralogous"
    }
    
    if ( input$colourVar == "AvgLength" ) { 
      df.sel$csel <- df.sel$AvgLength 
      colourtitle <- "Average length of locus"
    }
    
    
    # ifelse(input$categ),0.5/nlevels(Category),0
    
    p <- ggplot(df.sel) +
      geom_point( data=df.sel,
                  aes( x=xsel, y=ysel, colour = csel ), 
                  position = position_jitter(width = ifelse(input$categ,0.3,0)),
                  size=5, alpha=.5 ) +
      geom_hline( yintercept = mean(df.sel$ysel),
                  size=1, colour="blue", linetype="dashed" ) +
      ggtitle("Genetic diversity of loci") + 
      xlab("") +
      coord_flip() +
      theme_minimal() + 
      theme( axis.text.y = element_text(size=ifelse(input$categ,14,0)),
             axis.ticks  = element_line(size=ifelse(input$categ,0.5,0)),
             plot.title = element_text(face="bold"),
             axis.title.x = element_text(vjust=-.5, size=14),
             legend.position = "bottom", legend.key.width = unit(50, "pt") ) 

      

    if ( is.numeric(df.sel$csel) ) { 
      p <- p + guides(size = FALSE, alpha = FALSE, colour = guide_colorbar(title.position = "top") ) +  #+ #colour=FALSE, 
        scale_color_gradient2(low="yellow", 
                                     mid="red", 
                                     high="blue", 
                                     midpoint = median(df.sel$csel), 
                                     guide_legend(title = colourtitle))
    }
    else {
      p <- p + scale_colour_discrete(guide_legend(title = colourtitle))
    }
    
    # adding labels for tagged points 
    if ( labeltag$use )
    {
      p <- p +
      geom_text(  data = subset (df.sel, Locus %in% labeltag$list),
                  size=4, alpha=.8, vjust=-.5, angle = 30,
                  aes( x = xsel, y = ysel,
                        label= as.character(Locus) ) )
    }
    
    
    # putting in correct yVar label depending on variable chosen 
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
  
  
  # what actually does the plot, which also saves a PDF of it in case download is needed 
  output$mainplot <- renderPlot({
    return ( scatterPlot() )
  })
  
  # prints the download button and knows what file to make available 
  output$downloadPlot <- downloadHandler(
    filename = function() 
      { "plot.pdf" },
    content = function(file)
      { ggsave(file, width = 8, height = 8, plot = scatterPlot(), device = "pdf") }
  )
  
})