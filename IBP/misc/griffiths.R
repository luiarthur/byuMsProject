# Algorithm:
#   Start with arbitrary binary matrix:
#   For each row i,
#     For each column k,
#
#       if m-ik > 0, set zik =0 with prob eq.18
#       else delete that column
#
#       At the end of the row, add Poisson(a/N) new columns that have 1's
#     End 
#   End  
#
#  After sufficiently many passes through the rows, 
#  the resulting matrix will be a DRAW from the 
#  distribution P(Z) (eq.15)

rm(list=ls())
options("width"=120)

library(foreach)
library(doMC)
registerDoMC(16)

plot.binary <- function(M) {
  
  N <- nrow(M)
  K <- ncol(M)

  x <- rep(c(1:K),N)
  y <- rep(c(1:N),each=K)
  plot(x,y,xlim=c(1,K),ylim=rev(range(y)),
       ylab='Observation Number',xlab='Feature Number',
       main=paste('Feature Matrix'),pch=15,col='pink')

  colorColumns <- function(r){  
    p <- NULL
    for (c in 1:K){
      if (M[r,c]==1) {
        p <- rbind(c(c,r),p)
      } 
    }
    return(p)
  }

  pts <- foreach(n=1:N,.combine=rbind) %dopar% colorColumns(n)
  points(pts,pch=15,col='blue')
}

# Start with arbitrary binary matrix:
  N <- 50
  K <- 50
  Z <- matrix(sample(0:1,N*K,replace=T),nrow=N)
  plot.binary(Z)

  for (i in 1:N){
    k <- 1
    end.of.row <- F

    k <- 1
    while (!end.of.row){
       
    }
  }


