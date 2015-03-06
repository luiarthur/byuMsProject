source("ddibp.R")
library(xtable)

B <- 1e4 # problems
a <- .5

D1 <- matrix(c(0,1,9,
               1,0,1,
               9,1,0),3,3)

D2 <- matrix(c(0,1,1,
               1,0,9,
               1,9,0),3,3)

D3 <- matrix(c(0,1,1,
               1,0,9,
               1,9,0),3,3)

D4 <- matrix(c(0,1,9,
               1,0,1,
               9,1,0),3,3)

D5 <- matrix(c(0,9,1,
               9,0,1,
               1,1,0),3,3)

D6 <- matrix(c(0,9,1,
               9,0,1,
               1,1,0),3,3)

D71 <- matrix(c(0,1,9,
                1,0,9,
                9,9,0),3,3)

D7 <- matrix(c(0,9,1,
               9,0,9,
               1,9,0),3,3)

D8 <- matrix(c(0,9,9,
               9,0,1,
               9,1,0),3,3)

D9 <- matrix(c(0,1,1,9,
               1,0,1,1,
               1,1,0,1,
               9,1,1,0),4,4)

D10 <- matrix(c(0,1,1,9,
                1,0,9,1,
                1,9,0,1,
                9,1,1,0),4,4)

D11 <- matrix(c(0,1,1,9,
                1,0,9,1,
                1,9,0,1,
                9,1,1,0),4,4)

D12 <- matrix(c(0,9,1,9,
                9,0,9,9,
                1,9,0,9,
                9,9,9,0),4,4)

D13 <- matrix(c(0,9,9,1,
                9,0,9,9,
                9,9,0,9,
                1,9,9,0),4,4)

D99 <- matrix(9,7,7)
diag(D99) <- 0
D99[1,7] <- D99[7,1] <- 1

D55 <- matrix(9,5,5)
diag(D55) <- 0
D55[1,5] <- D55[5,1] <- 1

D54 <- matrix(9,5,5)
diag(D54) <- 0
D54[1,4] <- D54[4,1] <- 1

D5666 <- matrix(9,5,5)
diag(D5666) <- 0
D5666[1,3] <- D5666[3,1] <- 1
D5666[2,5] <- D5666[5,2] <- 1


#D <- matrix(0,3,3); D[which(lower.tri(D))] <- 1; 
#D <- matrix(1,3,3); D[which(upper.tri(D))] <- 0; 
#D <- matrix(c(0,1,20,1000,
#              1,0,50,30,
#              20,50,0,4,
#              1000,30,4,0),4,4)

#TEST: Permutations?
#Z <- lapply(as.list(1:B),function(x) {ot <- Sys.time()
#                                      o <- raibp(N=nrow(D55),a=2,D=D55,l=exp.decay,
#                                                 perm=T)
#                                      count.down(ot,x,B); o})
#EZ <- sum.matrices(Z)/B

exp.decay <- function(s,t,d) ifelse(s>t,0,exp(-d[s,t]))
exp.f <- function(d) exp(-d)
inv <- function(s,t,d) ifelse(s>t,0,1/d[s,t])
inv.f <- function(d) 1/d

one.sim <- function(D.name,a,B=1e4,num.cex=1,printProgress=F,lF=function(x) 1) {
  D <- eval(parse(text=D.name))

  cat("Getting Draws (1/4): \n")
  Z.o <- lapply(as.list(1:B),function(x) {ot <- Sys.time() 
                                          o <- raibp(N=nrow(D),a=a,lf=lF)
                                          count.down(ot,x,B); o})
  Z.d <- lapply(as.list(1:B),function(x) {ot <- Sys.time(); o <- rddibp(a=a,D=D)
                                          count.down(ot,x,B); o}) 
  Z.a <- lapply(as.list(1:B),function(x) {ot <- Sys.time()
                                          o <- raibp(N=nrow(D),a=a,D=D,l=exp.decay,
                                                     lf=lF)
                                          count.down(ot,x,B); o})
  #Z.a <- lapply(as.list(1:B),function(x) raibp(N=nrow(D),a=a,D=D,l=inv.decay))
  #Z.a <- lapply(as.list(1:B),function(x) raibp(N=nrow(D),a=a))

  cat("Calculating Expected Values (2/4): \n")
  EZO <- sum.matrices(Z.o)/B
  EZD <- sum.matrices(Z.d)/B
  EZA <- sum.matrices(Z.a)/B

  # Expected number of columns are different for ddIBP.
  mncolo <- mean(unlist(lapply(Z.o,ncol)))
  mncold <- mean(unlist(lapply(Z.d,ncol)))
  mncola <- mean(unlist(lapply(Z.a,ncol)))


  cat("Getting Unique Matrices from Draws (3/4): \n")
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

    cmo <- get.freqs(Z.o)
    cma <- get.freqs(Z.a)
    cmd <- get.freqs(Z.d)


    if (printProgress) cat("\tCalculating Info[,1/7]\n")
    info[,2] <- unlist(lapply(u.all$m,function(z) get.freq(z,cmo)))
    if (printProgress) cat("\tCalculating Info[,2/7]\n")
    info[,3] <- unlist(lapply(u.all$m,function(z) daibp(z,a=a)))
    if (printProgress) cat("\tCalculating Info[,3/7]\n")
    info[,4] <- unlist(lapply(u.all$m,function(z) get.freq(z,cma)))
    if (printProgress) cat("\tCalculating Info[,4/7]\n")
    info[,5] <- unlist(lapply(u.all$m,function(z) daibp(z,a=a,D=D,l=exp.decay)))
    if (printProgress) cat("\tCalculating Info[,5/7]\n")
    info[,6] <- unlist(lapply(u.all$m,function(z) get.freq(z,cmd)))
    info[,7] <- NA
    if (printProgress) cat("\tCalculating Info[,6/7]\n")
    info <- info[order(info[,7],decreasing=TRUE),]
    info[,1] <- 1:n
    if (printProgress) cat("Done Calculating Info\n")

    list("info"=info,"unique"=u.all)
  }

  cat("Comparing Methods (4/4): \n")
  M <- compare()
  cat("Done\n")

  options("width"=80)
  list("EZO"=EZO,"EZD"=EZD,"EZA"=EZA,"uzo"=uzo,"uzd"=uzd,"uza"=uza,"M"=M,
       "mncolo"=mncolo,"mncola"=mncola,"mncold"=mncold,"D"=D,"a"=a,"B"=B,
       "Zo"=Z.o,"Za"=Z.a,"Zd"=Z.d)
}

##Comment out:
## Graphs:
#source("ddibp.R")
#a <- 1.5
#result <- one.sim("D54",a=a,B=10000,printProg=T,lF=exp.f)
#result <- one.sim("D55",a=a,B=10000,printProg=T,lF=exp.f)
#result <- one.sim("D71",a=a,B=10000,printProg=T,lF=exp.f)
#result <- one.sim("D7",a=a,B=10000,printProg=T,lF=exp.f)
#result <- one.sim("D8",a=a,B=10000,printProg=T,lF=exp.f)
#result <- one.sim("D5666",a=2,B=10000,printProg=T,lF=exp.f)
#
###X11()
#pdf("../../../prospectus/images/eSim.pdf")
#  par(mfrow=c(3,1))
#    a.image(result$EZO,number=T,main=paste("E[IBP], E[ncol] =",result$mncolo),
#            num.cex=1)
#    a.image(result$EZA,number=T,main=paste("E[AIBP], E[ncol] =",result$mncola),
#            num.cex=1)
#    a.image(result$EZD,number=T, main=paste("E[ddIBP], E[ncol] =",
#            result$mncold),num.cex=.7)
#  par(mfrow=c(1,1))
#dev.off()
#
#sink("../../../prospectus/images/D54.tex")
#  print(xtable(D54,digits=0,align=rep("",ncol(D54)+1)),include.rowname=F,
#        include.colname=F,hline.after=F,hline=F,tabular.environment="pmatrix")
#sink()

#B <- 5e4
#X <- matrix(c(1,1,0,0,0,0,
#              0,0,1,1,0,0,
#              1,1,0,0,0,0,
#              0,0,0,0,1,1,
#              0,0,0,0,0,0),5,6,byrow=T)
#
#Zs <- lapply(as.list(1:B),function(x) {
#        ot <- Sys.time()
#      # o <- F.(matrix(X[1,1:2],nrow=1),a=2,start=2,end=5,
#      #         lam=function(s,t,d=D5666) exp(-d[s,t]))
#        o <- F.(X[1:2,1:4],a=2,start=3,end=5,lam=function(s,t,d=D5666) exp(-d[s,t]))
#      # o <- F.(X[1:3,1:4],a=2,start=4,end=5,lam=function(s,t,d=D5666) exp(-d[s,t]))
#      # o <- F.(X[1:4,1:6],a=2,start=5,end=5,lam=function(s,t,d=D5666) exp(-d[s,t]))
#        count.down(ot,x,B); o
#      })
#
#EZs <- sum.matrices(Zs)/B

########################
#B <- 5e4
#result<- one.sim("D5666",B=B,a=2,printProgress=T)
#Z.a <- result$Za
#Z.o <- result$Zo
#Z.d <- result$Zd
#
#calc.settings <- function(Zs) {
#  Z11 <- lapply(Zs,function(z) if (sum(z[1,])==2) z)
#  Z11 <- Z11[!unlist(lapply(Z11, is.null))]
#  EZ11 <- sum.matrices(Z11)/length(Z11)
#  u.z11 <- unique.matrix(Z11)
#
#  Z11.0011 <- lapply(Z11,function(z) if (all(z[2,1:2]==c(0,0)) && sum(z[2,])==2) z)
#  Z11.0011 <- Z11.0011[!unlist(lapply(Z11.0011,is.null))]
#  EZ11.0011 <- sum.matrices(Z11.0011)/length(Z11.0011)
#  print(length(Z11.0011))
#
#  EZ11.0011
#}
#
#Eo <- calc.settings(Z.o) # 630
#Ea <- calc.settings(Z.a) # 584
#Ed <- calc.settings(Z.d) # 3678
#
#
#pdf("../../../prospectus/images/eSimFixed.pdf")
#  par(mfrow=c(3,1))
#    a.image(round(Eo,5),numbers=T,num.cex=.8,main="E [IBP(2) | First 2 Rows]")
#    a.image(round(Ea,5),numbers=T,num.cex=.8,main="E [AIBP(2) | First 2 Rows,D]")
#    a.image(round(Ed,5),numbers=T,num.cex=.7,main="E [ddIBP(2) | First 2 Rows,D]")
#  par(mfrow=c(1,1))
#dev.off()

#uni <- one.sim("Duni",B=B,a=.5,printProgress=T)
#temp <- one.sim("Duni",B=B,a=.5,printProgress=T)
#pdf("../../../prospectus/images/uni.pdf")
#  par(mfrow=c(3,1))
#    a.image(round(uni$EZO,5),numbers=T,num.cex=.8,main="E [IBP(.5)]")
#    a.image(round(uni$EZA,5),numbers=T,num.cex=.8,main="E [AIBP(.5) | f=1]")
#    a.image(round(temp$EZO,5),numbers=T,num.cex=.8,main="E [ddIBP(.5) | f=1]")
#  par(mfrow=c(1,1))
#dev.off()

