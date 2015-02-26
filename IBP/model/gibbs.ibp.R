pb <- function(i,n) { 
  cat(paste0("\rProgress: ",round(i*1000/n)/10,"%"))
  if (i==n) cat("\n")
}

gibbs.ibp <- function(a=1,N=3,B=1e3,burn=B*.1) {
  
  Z <- list()
  length(Z) <- B
  Z[[1]] <- matrix(0,N,0)
  #Z[[1]] <- matrix(1,N,20)
  
  P <- function(z,i,k) { # P(z_{ik}=1 | z_{-i,k})
    sum(z[-i,k]) / N
  }

  for (b in 2:B) { # number of gibbs iterations
    pb(b,B)
    z <- Z[[b-1]]

    for (i in 1:N) { # iterate through all rows of Z
      K <- ncol(z)
      
      k <- 1
      while (k <= K) { # iterate through all columns of Z
        if (K>0) {
          p <- P(z,i,k)   
          u <- runif(1)
          if (p==0) {
            z <- as.matrix(z[,-k])
            k <- k-1
          } else if (p>u){
            z[i,k] <- 1
          } else {
            z[i,k] <- 0
          }
        }

        k <- k + 1
        K <- ncol(z)
      }

      new <- rpois(1,a/N)
     
      if (new>0) {
        col0 <- matrix(0,N,new)
        col0[i,] <- 1
        z <- as.matrix(cbind(z,col0))
      }

    }
    
    Z[[b]] <- z
  }

  out <- Z[(burn+1):B]
  out
}


## To check if my sampler is correct, I make sure Na, the expected number of counts
## equals N * alpha (or N * a)
#a <- .01
#N <- 3
#B <- 1e5
#Z <- gibbs.ibp(a=a,N=N,B=B)
#samp.counts <- unlist(lapply(Z,sum))
#(Na <- mean(samp.counts)) # Na = expected counts in Z
#plot(samp.counts,type="l",ylim=c(min(samp.counts),max(samp.counts)*1.1),
#     main=paste0("Gibbs: Trace Plot of Sample Counts, a=",a,", N=",N,", B=",B))
#abline(h=Na,col="orange",lwd=4,lty=1)
#abline(h=N*a,col="blue",lwd=4,lty=3)
#legend("topright",legend=c(paste0("Mean of Empirical Counts=",round(Na,4)),
#                           paste0("Theoretical Mean Na=",N*a)),
#                           col=c("orange","blue"),lwd=4,lty=c(1,3))

##source("model.R") # dIBP
#source("expectations.R") # rIBP
#source("sumMatrices.R")
##source("../IBP/IBP.R")
#
#
#library(doMC); registerDoMC(strtoi(system("nproc",intern=T))/2)
#L <- foreach(i=1:B) %dopar% rIBP(a=a,ppl=N)
#plot(unlist(lapply(L,sum)),type="l")
#mean(unlist(lapply(L,sum)))
#
#um.Z <- unique.matrix(Z)
#um.L <- unique.matrix(L)
#
#um.Z[[1]] /(B*.9)
#um.L[[1]] / B

