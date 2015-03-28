source("aibp2.R")
source("countdown.R")
library(doMC) 
registerDoMC(system("nproc",intern=TRUE))

calc.settings <- function(Zs) {
  Z11 <- lapply(Zs,function(z) if (sum(z[1,])==2) z)
  Z11 <- Z11[!unlist(lapply(Z11, is.null))]
  EZ11 <- sum.matrices(Z11)/length(Z11)
  u.z11 <- unique.matrix(Z11)

  Z11.0011 <- lapply(Z11,function(z) if (all(z[2,1:2]==c(0,0)) && sum(z[2,])==2) z)
  Z11.0011 <- Z11.0011[!unlist(lapply(Z11.0011,is.null))]
  EZ11.0011 <- sum.matrices(Z11.0011)/length(Z11.0011)
  print(length(Z11.0011))

  EZ11.0011
}

B <- 1e5
D <- matrix(9,5,5)
diag(D) <- 0
D[1,3] <- D[3,1] <- D[2,5] <- D[5,2] <- 1

D2 <- matrix(9,5,5)
diag(D2) <- 0
D2[1,3] <- D2[3,1] <- D2[2,5] <- D2[5,2] <- D2[3,4] <- D2[4,3] <- 1 

D6 <- matrix(9,6,6)
diag(D6) <- 0
D6[1,3] <- D6[3,1] <- D6[2,6] <- D6[6,2] <- D6[3,4] <- D6[4,3] <- 1


Zo <- foreach(i=1:B) %dopar% {ot <- Sys.time()
                              o <- raibp(a=2,N=6,l=function(s,t,d) exp(-d[s,t]))
                              count.down(ot,i,B); o}
Zop <- foreach(i=1:B) %dopar% {ot <- Sys.time()
                               o <- raibp(a=2,N=6,l=function(s,t,d) 
                                          exp(-d[s,t]),perm=TRUE)
                               count.down(ot,i,B); o}
Zs <- foreach(i=1:B) %dopar% {ot <- Sys.time()
                              o <- raibp(a=2,D=D6,l=function(s,t,d) exp(-d[s,t]))
                              count.down(ot,i,B); o}
Zp <- foreach(i=1:B) %dopar% {ot <- Sys.time()
                              o <- raibp(a=2,D=D6,l=function(s,t,d) 
                                         exp(-d[s,t]),perm=TRUE)
                              count.down(ot,i,B); o}           

EO <- sum.matrices(Zo) / B
EOp <- sum.matrices(Zop) / B
EZ <- sum.matrices(Zs) / B
EP <- sum.matrices(Zp) / B

EOc <- calc.settings(Zo) #128
EZc <- calc.settings(Zs) #128
EPc <- calc.settings(Zp) #165

par(mfrow=c(3,2))
  a.image(round(EO,3),numbers=TRUE,main="IBP")
  a.image(round(EOc,3),numbers=TRUE,main="IBP")
  a.image(round(EZ,3),numbers=TRUE,main="Not Permuted")
  a.image(round(EZc,3),numbers=TRUE,main="Not Permuted")
  a.image(round(EP,3),numbers=TRUE,main="Permuted")
  a.image(round(EPc,3),numbers=TRUE,main="Permuted")
par(mfrow=c(1,1))

# THIS IS ALSO IMPORTANT: Make sure the permuted version reduces to the IBP
# when the observations are equidistant.
#par(mfrow=c(2,1))
#  a.image(round(EO,3),numbers=TRUE,main="IBP Not Permuted")
#  a.image(round(EOp,3),numbers=TRUE,main="IBP Permuted")
#par(mfrow=c(1,1))

up <- unique.matrix(Zp)
which(up[[1]] > B*.001)
up[[1]][1:10]
up[[2]][1:10]
