library(shiny)

shinyUI(fluidPage(
  titlePanel("GbGDiv Stats & Graphs"),

  sidebarLayout(

    sidebarPanel(
      
      h4("Upload your table of GbGDiv results"),
      
      fileInput('resultstable', 'Choose ResultsTable file',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      fileInput('categorytable', 'Choose Category file',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      conditionalPanel(condition = "output.catfileUp == true",
                       checkboxInput(inputId = "categ",
                                     label = "Use categories",
                                     value = FALSE)),
                       
      tags$hr(),
      
      h4("Choose data set"),
      
      selectInput(inputId = "yVar",
                  label = "Which varible do you want to explore?",
                  choices = c("AllelicDiv","ADivNM","VSitesNuc","VSitesAA","RatioCount","RatioVS"),
                  selected = "AllelicDiv"),      
      
      tags$hr(),
      
      h4("Options"),
        
      h6("Label axes with standard deviations, not absolute values"),
      checkboxInput("zscore", label = "Use z-scores", value = FALSE)
      

    ), # closes sidebar

    mainPanel(
      plotOutput(outputId = "mainplot" )
      #tableOutput("my_output_data")
      
    ) # closes main
  ) # closes layout
)) #closes fluidPage & shinyUI