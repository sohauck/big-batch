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
      
      sliderInput("percexc", "Minimum percentage of isolates tagged for locus to be included:", 
                  min = 0, max = 100, value = 10, step = 1),

      tags$hr(),
      
      h4("Choose variable to plot"),
      
      selectInput(inputId = "yVar",
                  label = "Which varible do you want to explore?",
                  choices = c("AllelicDiv","ADivNM","VSitesNuc","VSitesAA","RatioCount","RatioVS"),
                  selected = "AllelicDiv"),      
      
      tags$hr(),
      
      h4("Optional parameters"),
        
      checkboxInput("zscore", label = "Use z-scores for y-axis values", value = FALSE),
      
      fileInput('categorytable', 'Choose Category file',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      conditionalPanel(condition = "output.catfileUp == true",
                       checkboxInput(inputId = "categ",
                                     label = "Use categories",
                                     value = FALSE)),
      
      tags$hr(),
      
      h4("Download images")
      
      # output$downloadData <- downloadHandler(
      #   filename = function() { 
      #     paste(input$yVar, '.csv', sep='') 
      #   },
      #   content = function(file) {
      #     ( ggsave("GbGDiv.png", height = 9, width = 12, dpi = 100), file )
      #   }
      
    ), # closes sidebar

    mainPanel(
      tabsetPanel(type = "tabs", 
        tabPanel("Scatter Plot", plotOutput(outputId = "mainplot")), 
        tabPanel("Distribution Plot", plotOutput(outputId = "distplot"), plotOutput(outputId = "corrplot")),
        tabPanel("Data Table", tableOutput("data.table"))
       )

    ) # closes main
  ) # closes layout
)) #closes fluidPage & shinyUI