library(shiny)

shinyUI(fluidPage(
  tags$head(includeScript("google-analytics.js")),
  
  titlePanel("GbGDiv Stats & Graphs + (Google Analytics!)"),

  sidebarLayout(

    sidebarPanel(
      
      h4("Upload your GbGDiv table"),
      
      fileInput('resultstable', 'Choose ResultsTable file',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      sliderInput("percexc", "Max 'missing' allele tag percentage per locus:", 
                  min = 0, max = 100, value = 5, step = .1),

     # textOutput("exclusions"), 
      tagAppendAttributes(textOutput("exclusions"), style="white-space:pre-wrap;"),
      
      tags$hr(),
      
      h4("Choose variables to plot"),
      
      selectInput(inputId = "yVar",
                  label = "Which variable do you want to explore?",
                  choices = c("AllelicDiv","ADivNM","VSitesNuc","VSitesAA","RatioCount","RatioVS"),
                  selected = "AllelicDiv"),      
      
     selectInput(inputId = "colourVar",
                 label = "Which variable do you want in the point colours?",
                 choices = c("None (labels only)","Missing","Paralogous","AvgLength"),
                 selected = "None"),      
     
      tags$hr(),
      
      h4("Optional parameters"),
        
      checkboxInput("zscore", label = "Use z-scores for y-axis values", value = FALSE),
      
     
      conditionalPanel(condition = "output.catinc == true",
                       checkboxInput(inputId = "categ",
                                     label = "Use categories",
                                     value = FALSE)),
      
      actionButton('resetlabel', 'Reset labelled points'),
      
      tags$hr(),
      
      # download from Scatter Plot
      conditionalPanel(condition = "input.conditionaltab == 1",
                       downloadButton('downloadPlot', 'Download Plot')
                       )

      
    ), # closes sidebar

    mainPanel(
      tabsetPanel(id = "conditionaltab", type = "tabs", 
        tabPanel("Scatter Plot", value = "1", 
                 plotOutput(outputId = "mainplot", height = 700,
                            dblclick = "mp_dblclick",
                            brush = brushOpts(
                              id = "mp_brush"))
                 #verbatimTextOutput("brushinfo")
                 ), 
        tabPanel("Excluding points", value = "2",
                 helpText("Choose what percentage of isolates a locus can be marked 'missing' (due to no allele designation)",
                          "in and still included in the analysis. In other words, a max 'missing' percentage.",
                          "Use the slider on the left to choose a cutoff that excludes loci that have been",
                          "measured as having low diversity (low y-values in the bottom plot) when they have",
                          "high 'missing' counts (high x-values in the same). The upper plot shows the distribution",
                          "of loci over the range of 0 to the total count of isolates, where the red-dotted line indicated that all",
                          "loci to its right have been excluded from all other graphs and tables in this application."),
                 plotOutput(outputId = "distplot", height = 400), 
                 plotOutput(outputId = "corrplot")),
        tabPanel("Data Table", value = "3", 
                 dataTableOutput("data.table"))
       )

    ) # closes main
  ) # closes layout
)) #closes fluidPage & shinyUI