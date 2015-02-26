source("../ddibp.R",chdir=T)
library(xtable)
library(shiny)

B <- 1e5 # problems
a <- .5

# 1 close to 2, 2 close-ish to 3, 1 far from 3. 
D1 <- matrix(c(0,1,10,
               1,0,3,
               10,3,0),3,3)

# 1 far from to 2, 2 close-ish to 3, 1 close from 3. 
D2 <- matrix(c(0,10,1,
               10,0,3,
               1,3,0),3,3)

# 1 close to 2, 2 far from 3, 1 close-ish 3. 
D3 <- matrix(c(0,1,3,
               1,0,10,
               3,10,0),3,3)

#D <- matrix(0,3,3); D[which(lower.tri(D))] <- 1; 
#D <- matrix(1,3,3); D[which(upper.tri(D))] <- 0; 
#D <- matrix(c(0,1,20,1000,
#              1,0,50,30,
#              20,50,0,4,
#              1000,30,4,0),4,4)

exp.decay <- function(s,t,d) ifelse(s>t,0,exp(-d[s,t]))
inv <- function(s,t,d) ifelse(s>t,0,1/d[s,t])

one.sim <- function(D,D.name,a,B=1e4,num.cex=1) {
  cat("Getting Draws (1/6): \n")
  Z.o <- lapply(as.list(1:B),function(x) {ot <- Sys.time() 
                                          o <- raibp(N=nrow(D),a=a)
                                          count.down(ot,x,B); o})
  Z.d <- lapply(as.list(1:B),function(x) {ot <- Sys.time(); o <- rddibp(a=a,D=D)
                                          count.down(ot,x,B); o}) 
  Z.a <- lapply(as.list(1:B),function(x) {ot <- Sys.time()
                                          o <- raibp(N=nrow(D),a=a,D=D,l=exp.decay)
                                          count.down(ot,x,B); o})
  #Z.a <- lapply(as.list(1:B),function(x) raibp(N=nrow(D),a=a,D=D,l=inv.decay))
  #Z.a <- lapply(as.list(1:B),function(x) raibp(N=nrow(D),a=a))

  cat("Calculating Expected Values (2/6): \n")
  EZO <- sum.matrices(Z.o)/B
  EZD <- sum.matrices(Z.d)/B
  EZA <- sum.matrices(Z.a)/B

  # Expected number of columns are different for ddIBP.
  mncolo <- mean(unlist(lapply(Z.o,ncol)))
  mncold <- mean(unlist(lapply(Z.d,ncol)))
  mncola <- mean(unlist(lapply(Z.a,ncol)))

  cat("Plotting Graphs (3/6): \n")
  png(paste0("www/",D.name,".png"),height=700)
  #png(paste0("pdf/",D.name,".png"))
    par(mfrow=c(4,1))
      #nr <- nrow(D)
      f <- exp(-D)
      f[which(upper.tri(f))] <- 0
      #f <- matrix(prettyNum(f,width=11,format="fg"),nr,nr)
      a.image(round(f,5),main="Proximity Matrix",axis=F,numb=T,num.cex=num.cex)
      a.image(EZO,main=paste("IBP, E[ncol] =",mncolo),cex.axis=1,numb=T,
              num.cex=num.cex)
      a.image(EZD,main=paste("ddIBP, E[ncol] =",mncold),cex.axis=1,numb=T,
              num.cex=num.cex*.8)
      a.image(EZA,main=paste("aIBP, E[ncol] =",mncola),cex.axis=1,numb=T,
              num.cex=num.cex)
    par(mfrow=c(1,1))
  dev.off()

  cat("Getting Unique Matrices from Draws (4/6): \n")
  uzo <- unique.matrix(Z.o)
  uzd <- unique.matrix(Z.d)
  uza <- unique.matrix(Z.a)

  # make a table:
  compare <- function() {

    options("width"=180)
    u.all <- unique.matrix(c(Z.o,Z.d,Z.a))
    n <- length(u.all$c)
    info <- matrix(0,n,7)
    colnames(info) <- c("matrix#",paste(rep(c("IBP","aIBP","ddIBP"),each=2),
                        c("empirical","theoretical")))

    # Make this faster using unique.matrix
    #get.freq <- function(m,Zs) {
    #  count <- 0
    #  ll <- length(Zs)

    #  for (i in 1:ll) {
    #    z <- Zs[[i]]
    #    if (ncol(m)==ncol(z) && nrow(m)==nrow(z)) {
    #      if (ncol(m)==0 || all(m==z)){
    #        count <- count + 1
    #      }
    #    }
    #  }
    #  
    #  count/B
    #}
    cmo <- get.freqs(Z.o)
    cma <- get.freqs(Z.a)
    cmd <- get.freqs(Z.d)

    info[,1] <- 1:n
    info[,2] <- unlist(lapply(u.all$m,function(z) get.freq(z,cmo)))
    info[,3] <- unlist(lapply(u.all$m,function(z) daibp(z,a=a)))
    info[,4] <- unlist(lapply(u.all$m,function(z) get.freq(z,cma)))
    info[,5] <- unlist(lapply(u.all$m,function(z) daibp(z,a=a,D=D,l=exp.decay)))
    info[,6] <- unlist(lapply(u.all$m,function(z) get.freq(z,cmd)))
    info[,7] <- NA
    head(info,20)

    list("info"=info,"unique"=u.all)
  }

  cat("Comparing Methods (5/6): \n")
  M <- compare()

  cat("Printing output (6/6): \n")
  sink(paste0("www/out",D.name,".txt"))
    print(M)
  sink()

  options("width"=80)
  list("EZO"=EZO,"EZD"=EZD,"EZA"=EZA,"uzo"=uzo,"uzd"=uzd,"uza"=uza,"M"=M,
       "mncolo"=mncolo,"mncola"=mncola,"mncold"=mncold)
}


out.D1 <- one.sim(D1,"D1",a=.5,B=B,num.cex=1)
out.D2 <- one.sim(D2,"D2",a=.5,B=B,num.cex=1)
out.D3 <- one.sim(D3,"D3",a=.5,B=B,num.cex=1)


createTable <- function (D.name,num.cex=c(1,1,1,1)) {
  matfreq <- eval(parse(text=paste0("out.",D.name,"$M$info")))
  um <- eval(parse(text=paste0("out.",D.name,"$M$unique$matrix")))
  D <- eval(parse(text=D.name))

  mncolo <- eval(parse(text=paste0("out.",D.name,"$mncolo")))
  mncola <- eval(parse(text=paste0("out.",D.name,"$mncola")))
  mncold <- eval(parse(text=paste0("out.",D.name,"$mncold")))


  EZO <- eval(parse(text=paste0("out.",D.name,"$EZO")))
  EZA <- eval(parse(text=paste0("out.",D.name,"$EZA")))
  EZD <- eval(parse(text=paste0("out.",D.name,"$EZD")))
  

  png(paste0("www/",D.name,".png"),height=700)
  #png(paste0("pdf/",D.name,".png"))
    par(mfrow=c(4,1))
      #nr <- nrow(D)
      f <- exp(-D)
      f[which(upper.tri(f))] <- 0
      #f <- matrix(prettyNum(f,width=11,format="fg"),nr,nr)
      a.image(round(f,5),main="Proximity Matrix",axis=F,numb=T,num.cex=num.cex[1])
      a.image(EZO,main=paste("IBP, E[ncol] =",mncolo),cex.axis=1,numb=T,
              num.cex=num.cex[2])
      a.image(EZD,main=paste("ddIBP, E[ncol] =",mncold),cex.axis=1,numb=T,
              num.cex=num.cex[3])
      a.image(EZA,main=paste("aIBP, E[ncol] =",mncola),cex.axis=1,numb=T,
              num.cex=num.cex[4])
    par(mfrow=c(1,1))
  dev.off()


  runApp(list(
    ui = basicPage(
      h2("Matrix Frequencies"),
      sidebarLayout(position="right",
        sidebarPanel(
          #numericInput("distMat",label="Distance Matrix",value=1,min=1,max=)
          imageOutput("Image",height=700),
          br(),
          numericInput("matNum",label="Matrix #:",value=1,min=1,max=length(um)),
          uiOutput("matrix")
        ),
        mainPanel( 
          dataTableOutput("matFreq")
        )
      )
    ),
    server = function(input,output) {
      output$Image = renderImage({
        list(src=paste0("www/",D.name,".png"))
      }, deleteFile=FALSE)

      output$matFreq = renderDataTable({matfreq})

      output$matrix = renderUI({
        M <- um[[input$matNum]]
        if (ncol(M)>0) {
          M <- print(xtable(M, align=rep("", ncol(M)+1)),
                     floating=FALSE, tabular.environment="bmatrix", comment=FALSE, 
                     print.results=FALSE,hline.after=NULL,
                     include.rowname=FALSE,include.colname=FALSE)
          html <- paste0("$$", M, "$$")
          withMathJax(helpText(html))
        }
      })
    }
  ))
}

createTable("D1",num.cex=c(2,1.5,1.3,1.3))
createTable("D2")
createTable("D3")

# Want to click on matrix number and show matrix in line.
