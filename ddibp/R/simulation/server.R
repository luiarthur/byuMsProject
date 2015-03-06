printMatrix <- function(M) {
  if (ncol(M)>0) {
    M <- print(xtable(M, align=rep("", ncol(M)+1)),
               floating=FALSE, tabular.environment="bmatrix", comment=FALSE, 
               print.results=FALSE,hline.after=NULL,
               include.rowname=FALSE,include.colname=FALSE)
    html <- paste0("$$", M, "$$")
    withMathJax(helpText(html))
  } else {
    "<Empty Vector>"
  }
}

shinyServer(function(input,output) {
  source("ddibp.R",chdir=T)
  source("oneSim.R")
  library(xtable)
  
  #B <- reactive({input$its})
  #a <- reactive({input$alpha})
  #result <- NULL
  #submit <- input$submit
  
  result <- reactive({
    if (input$submit==0) return(NULL)
    one.sim(paste0("D",input$distMatNum),a=input$alpha,B=input$its,printProg=T)
  })

  output$distMat = renderUI({
    DM <- eval(parse(text=paste0("D",input$distMatNum))) # Distance Matrix
    D <- print(xtable(DM, align=rep("", ncol(DM)+1), digits=0),
               floating=FALSE, tabular.environment="pmatrix", comment=FALSE, 
               print.results=FALSE,hline.after=NULL,
               include.rowname=FALSE,include.colname=FALSE)
    html <- paste0("$$", D, "$$")
    withMathJax(helpText(html))
  })

  #output$Image = renderImage({
  #  list(src=paste0("www/D",input$distMatNum,".png"))
  #}, deleteFile=FALSE)

  output$matFreq = renderDataTable({
    if (input$submit==0) return(NULL)
    (result())$M$info
  })

  output$expVal = renderPlot({
    if(input$submit==0) return(NULL)
    par(mfrow=c(3,1))
      a.image(result()$EZO,number=T,main=paste("E[IBP], E[ncol] =",result()$mncolo))
      a.image(result()$EZA,number=T,main=paste("E[AIBP], E[ncol] =",result()$mncola))
      a.image(result()$EZD,number=T, main=paste("E[ddIBP], E[ncol] =",
              result()$mncold))
    par(mfrow=c(1,1))
  })

  output$matrix=renderUI({ 
    if(input$submit==0) return(NULL)
    printMatrix(result()$M$unique$matrix[[input$matNum]])
  })
})
