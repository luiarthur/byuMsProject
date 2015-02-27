shinyUI(fluidPage(
  titlePanel("Simulation Study: IBP vs. aIBP vs. ddIBP"),

  sidebarLayout(position="right",
    sidebarPanel(
      numericInput("its",label="Number of Iterations",value=1000),
      numericInput("alpha",label="Alpha",value=.5,min=1e-6),
      numericInput("matNum",label="Matrix #:",value=1,min=1),
      uiOutput("matrix"),
      numericInput("distMatNum",label="Distance Matrix",value=1,min=1,max=9),
      uiOutput("distMat"),
      tags$hr(), br(),
      
      tags$b("Expected Values:"),
      plotOutput("expVal")
    ),
    mainPanel( 
      dataTableOutput("matFreq")
    )
  )
))
