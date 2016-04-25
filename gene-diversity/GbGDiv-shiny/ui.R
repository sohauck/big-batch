library(shiny)

shinyUI(fluidPage(
  titlePanel("GbGDiv Stats & Graphs"),

  sidebarLayout(

    sidebarPanel(
      
      h4("Upload your GbGDiv table"),
      
      fileInput('resultstable', 'Choose ResultsTable file',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      sliderInput("percexc", "Percent of loci excluded due to high number of missing sequence tags:", 
                  min = 0, max = 50, value = 5, step = .1),

      textOutput("exclusions"), 
      
      tags$hr(),
      
      h4("Choose variable to plot"),
      
      selectInput(inputId = "yVar",
                  label = "Which varible do you want to explore?",
                  choices = c("AllelicDiv","ADivNM","VSitesNuc","VSitesAA","RatioCount","RatioVS"),
                  selected = "AllelicDiv"),      
      
      tags$hr(),
      
      h4("Optional parameters"),
        
      checkboxInput("zscore", label = "Use z-scores for y-axis values", value = FALSE),
      
      conditionalPanel(condition = "output.catinc == true",
                       checkboxInput(inputId = "categ",
                                     label = "Use categories",
                                     value = FALSE)),
      
      helpText("You can label any points by dragging over them and double clicking."),
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
                 plotOutput(outputId = "mainplot", height = 600,
                            dblclick = "mp_dblclick",
                            brush = brushOpts(
                              id = "mp_brush"))
                 #verbatimTextOutput("brushinfo")
                 ), 
        tabPanel("Excluding points", value = "2",
                 helpText("Choose what percentage of loci to exclude due to high counts of isolates",
                          "where the locus was marked 'missing' due to no allele designation.",
                          "Use the slider on the left to choose a cutoff that excludes loci that have been",
                          "measured as having low diversity (low y-values in the bottom plot) when they have",
                          "high 'missing' counts (high x-values in the same). The upper plot shows the distribution",
                          "of loci over the range of 'missing' values, where the red-dotted line indicated that all",
                          "loci to its right have been excluded from all other graphs and tables in this application."),
                 plotOutput(outputId = "distplot", height = 400), 
                 plotOutput(outputId = "corrplot")),
        tabPanel("Data Table", value = "3", 
                 dataTableOutput("data.table"))
       )

    ) # closes main
  ) # closes layout
)) #closes fluidPage & shinyUI