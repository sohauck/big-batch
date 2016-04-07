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
                  min = 0, max = 100, value = 10, step = 1),

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
      
      tags$hr(),
      
      #  && 
      #output.fileUp == true
      
      # download from Scatter Plot
      conditionalPanel(condition = "input.conditionaltab == 1",
                       downloadButton('downloadPlot', 'Download Plot'),
                       actionButton("exclude_toggle", "Toggle points"),
                       tags$br(),
                       actionButton("exclude_reset", "Reset")
      )

      
    ), # closes sidebar

    mainPanel(
      tabsetPanel(id = "conditionaltab", type = "tabs", 
        tabPanel("Scatter Plot", value = "1", 
                 plotOutput(outputId = "mainplot", 
                            click = "plot1_click",
                            brush = brushOpts(
                              id = "plot1_brush")),
                 verbatimTextOutput(outputId = "info")
                 ), 
        tabPanel("Distribution Plot", value = "2", 
                 plotOutput(outputId = "distplot", height = 250), 
                 plotOutput(outputId = "corrplot")),
        tabPanel("Data Table", value = "3", 
                 dataTableOutput("data.table"))
       )

    ) # closes main
  ) # closes layout
)) #closes fluidPage & shinyUI