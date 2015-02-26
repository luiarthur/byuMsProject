# FUNCTIONS: ########################################################
source("../IBP/model/sumMatrices.R")
source("../IBP/model/countdown.R")
library(doMC); registerDoMC(strtoi(system("nproc",intern=T))/2)
pb <- function(i,n) { 
  cat(paste0("\rProgress: ",round(i*1000/n)/10,"%"))
  if (i==n) cat("\n")
}

a.image <- function(Q,color=paste0("gray",100:0),...) {
  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
}
#####################################################################

a <- 3 # 12.RD
b <- 9 # 13.RL
c <- 1 # 23.DL

D <- matrix(c(0,a,b,
              a,0,c,
              b,c,0),3,3)

D1 <- matrix(c(0,1,1,
               1,0,1,
               1,1,0),3,3)


inv <- function(s,t,d=D) 1/d[s,t] # inverse distance metric

# Calculates Probability of Customer_i Getting Dish_k

f. <- function(x,i=2,lam=inv) {
  K <- ncol(x)
  f <- function(x,i,k) {
    out <- 0
    if (sum(x[,k]>0)) {
      ind <- which(x[,k]==1)
      out <- sum(lam(ind,i))
    }
    out
  }

  out <- 0
  if (K>0) {
    t <- apply(matrix(1:K),1,function(k) f(x,i,k)) 
    out <- rep(0,K)
    if (sum(t)>0) {
      out <- t/sum(t) * sum(x[1:(i-1),])/i
    }
  }
  out
}

rMyIBP <- function(N=3,a=3,D=NULL,l=inv) {
  K <- rpois(1,a)
  Z <- matrix(0,N,K) 
  Z[1,0:K] <- 1 # The first customer draws a POI(a) number of new dishes
  
  if (is.null(D)) {
    D <- matrix(1,N,N)
    diag(D) <- 0
  }
  L <- function(s,t,d) l(s,t,D)
 
  for (i in 2:N) {
    P <- f.(Z,i,lam=L)
    #print(P)
    if (K>0) Z[i,] <- P > runif(K)
    newK <- K+rpois(1,a/i)
    col0 <- matrix(0,N,newK-K)
    if (ncol(col0) > 0) {
      Z <- cbind(Z,col0)
      Z[i,(K+1):newK] <- 1
      K <- newK
    }
  }
  
  Z
}

# Simulate Properties:
B <- 100000
Zs <- foreach(i=1:B) %dopar% rMyIBP(a=3,D=D)
#Zs <- foreach(i=1:B) %dopar% rMyIBP(D=D1)
EZ <- sum.matrices(Zs) / length(Zs)
VZ <- sum.matrices(lapply(Zs,function(z)z^2))/length(Zs) - EZ^2
a.image(EZ) 
apply(EZ,1,sum) # should be a
sum(EZ) # should be aN
ncolZ <- unlist(lapply(Zs,ncol))
mean(ncolZ) # should be aHn = a*sum(1/1:N)

# Simulate X:
X <- matrix(c(1,1,1,0,
              0,1,0,1,
              0,0,0,0),3,4,byrow=T)

fx <- function(x=X,i=3,D,lam) {
  L <- function(s,t) lam(s,t,D)
  P <- f.(x,3,L)
  P > runif(4)
}

x3 <- foreach(i=1:B,.combine=rbind) %dopar% fx(X,3,D,inv)
pred <- apply(x3,2,mean)
pred.x3 <- X
pred.x3[3,] <- pred
a.image(pred.x3)
mean(apply(x3,1,sum))
pred

# What Next?
# Yes, my expected existing dishes = 5/3!
