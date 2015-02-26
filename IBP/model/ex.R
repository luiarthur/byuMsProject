ource("model.R")

f <- function(set){
  source("sets.R",local=T)
  print(n)
}

# 1) E[X|Z=z]
# Do n == d, n < d, and n > d

rxgz <- function(N=10000,set=1){
  source("sets.R",local=T)
  xgz <- foreach(i=1:N) %dopar% r.X.given.Z(Z,d,su,sv)
  list(xgz,a,Z)
}

xgz.gen <- function(set,N=1000,plotting=F){
  xgza <- rxgz(N,set=set)
  xgz <- xgza[[1]]; a <- xgza[[2]]; Z <- xgza[[3]]
  #(EXGZ <- Reduce('+', xgz) / N)
  #plot.expected.IBP(EXGZ,a,N=1)
  CXGZ <- Reduce('+', xgz)
  EXGZ <- CXGZ / N
  J <- matrix(1,ncol(Z),ncol(EXGZ))
  if (plotting){
    plot.expected.IBP(EXGZ,a,N=1,disNum=F,row.names=F,col.names=F,main="E[X|Z]")
  }
  # So, if Z[i,] has O_i 1's, then E[X|Z][i,] = Z%*%J[i,] = a vector of O_i's. 
  list("counts"=CXGZ,"EXGZ"=EXGZ,"Z"=Z,"J"=J)
}

#xgz.6 <- xgz.gen(6,10000)
#xgz.6$Z %*% xgz.6$J
#xgz.6$EXGZ
#
#xgz.7 <- xgz.gen(7,10000)
#xgz.7$Z %*% xgz.6$J
#xgz.7$EXGZ
#
#xgz.8 <- xgz.gen(8,10000)
#xgz.8$Z %*% xgz.6$J
#xgz.8$EXGZ
#
#xgz.9 <- xgz.gen(9,10000)
#xgz.9$Z %*% xgz.6$J
#xgz.9$EXGZ
#
#xgz.10 <- xgz.gen(10,10000)
#xgz.10$Z %*% xgz.6$J
#xgz.10$EXGZ



# 2) E[X] = E[E[X|Z]] # all cells tend to a
x.gen <- function(a,N=1000,plotting=F){
  n <- 3; d <- 2; su <- 1; sv <- 1
  x <- foreach(i=1:N) %dopar% r.X.given.Z.given.a(n,a,d,su,sv)

  CX <- Reduce('+', x) # all cells tend to a
  EX <- CX / N
  #J <- matrix(1,nrow(Z),ncol(EX))
  if (plotting) plot.expected.IBP(EX,5,N=N,row.n=F,col.n=F,main="E[X]=E[E[X|Z]]")

  list("counts"=CX,"EX"=EX)
  # Checking: (Correct): E[X] = E[E[X|Z]] = E[ZJ] = E[Z]J
  #EZ <- EIBP(a,N=N) / N
  #J <- matrix(1,ncol(EZ),d)
  #EZ%*%J
}

EZ <- function(n,a,d,N=1000) {
  ez <- expected.IBP(n=3,a=a,N=N) / N
  J <- matrix(1,ncol(ez),d)
  list("EZ"=ez,"J"=J,"EZJ"=ez%*%J)
}
#N <- 10000; a <- 3
#x.1 <- x.gen(3,N)
#EZJ <- EZ(n=3,a=3,d=ncol(x.1$EX),10000)
#EZJ$EZJ
#x.1$EX


