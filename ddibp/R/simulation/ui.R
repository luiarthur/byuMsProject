shinyUI(fluidPage(
  titlePanel("Simulation Study: IBP vs. aIBP vs. ddIBP"),

  sidebarLayout(position="right",
    sidebarPanel(
      numericInput("matNum",label="Matrix #:",value=1,min=1),
      uiOutput("matrix"),
      tags$hr(), br(),

      numericInput("its",label="Number of Iterations",value=1000),
      numericInput("alpha",label="Alpha",value=.5,min=1e-6),
      numericInput("distMatNum",label="Distance Matrix",value=1,min=1,max=9),
      uiOutput("distMat"),
      actionButton("submit","Submit"),
      tags$hr(), br(),

      tags$b("Expected Values:"),
      plotOutput("expVal")
    ),
    mainPanel( 
      dataTableOutput("matFreq")
    )
  )
))
