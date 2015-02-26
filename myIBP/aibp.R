# FUNCTIONS: ########################################################
source("../IBP/model/sumMatrices.R")
source("../IBP/model/countdown.R")
source("uniqueMatrix.R")
library(doMC); registerDoMC(strtoi(system("nproc",intern=T))/2)
pb <- function(i,n) { 
  cat(paste0("\rProgress: ",round(i*1000/n)/10,"%"))
  if (i==n) cat("\n")
}

a.image <- function(Q,color=paste0("gray",100:0),...) {
  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
}
#####################################################################
get.new.dish <- function(z) {
  N <- nrow(z)
  K <- ncol(z)
  x <- rep(0,N)
  x[1] <- sum(z[1,])
  for (i in 2:N) {
    if (sum(x[1:(i-1)])+1 <= K) x[i] <- sum(z[i,(sum(x[1:(i-1)])+1):K])
  }
  x
}

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
  if (is.null(K)) K <- 0
  f <- function(x,i,k) {
    out <- 0
    if (sum(x[1:(i-1),k]>0)) {
      ind <- which(x[1:(i-1),k]==1)
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


# For a GIVEN PERMUTATION!!!

raibp <- function(N=3,a=3,D=NULL,l=inv) {
  K <- rpois(1,a)
  Z <- matrix(0,N,K) 
  Z[1,0:K] <- 1 # The first customer draws a POI(a) number of new dishes
  
  # If no distance matrix is provided, customers will be equidistant.
  if (is.null(D)) {
    D <- matrix(1,N,N)
    diag(D) <- 0
  }

  if (N>=2) {
    for (i in 2:N) {
      P <- f.(Z,i,lam=function(s,t,d) l(s,t,D))
      #if (N>2) cat("PROB: ",paste(P),"\n")
      if (K>0) Z[i,] <- P > runif(K)

      newK <- K+rpois(1,a/i)
      col0 <- matrix(0,N,newK-K)

      if (ncol(col0) > 0) {
        Z <- cbind(Z,col0)
        Z[i,(K+1):newK] <- 1
        K <- newK
      }
    }
  }
  
  Z
}

ribp <- function(N=3,a=1){
  K <- rpois(1,a)
  Z <- matrix(0,N,K)
  Z[1,0:K] <- 1
  
  for (i in 2:N) {
    # Sample previously sampled dishes
    if (K>0) {
      mk <- apply(Z,2,sum)
      pk <- mk / i
      Z[i,] <- pk > runif(K)
    } 

    # Adding new dishes for customer i
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

dibp <- function(Z,a) {
  N <- nrow(Z)
  K <- ncol(Z)
  Hn <- sum(1/1:N)
  x <- get.new.dish(Z)
  mk <- apply(Z,2,sum)
  
  p <- a^K / prod(gamma(x+1)) * exp(-a*Hn) *
       prod(gamma(N-mk+1)*gamma(mk)/gamma(N+1))
  p
}

dibp2 <- function(Z,a) {
  N <- nrow(Z)
  K <- ncol(Z)
  Hn <- sum(1/1:N)
  x <- get.new.dish(Z)
  p <- a^K / prod(gamma(x+1)) * exp(-a*Hn)

  if (K>0) {
    for (i in 2:N) {
      p <- p / i^x[i]
      y <- sum(x[1:(i-1)])
      if (y>0) {
        for (k in 1:y) {# For each existing dish previous to customer i
          m_ik <- sum(Z[1:(i-1),k])
          p <- p * m_ik^Z[i,k] * (i-m_ik)^(1-Z[i,k]) / i
        }
      }
    }
  }

  p
}

daibp <- function(Z,a=3,D=NULL,l=inv,log=F) {
  N <- nrow(Z)
  K <- ncol(Z)
  x <- get.new.dish(Z)
  
  if (is.null(D)) {
    D <- matrix(1,N,N)
    diag(D) <- 0
  }

  p <- ifelse(log,0,1)
  if (K>0) {
    for (i in 1:N) {
      g <- dpois(x[i],a/i,log=log)
      h <- ifelse(log,0,1)
      existing.dishes <- sum(x[1:(i-1)])
      if (i>1 && existing.dishes>0) {
        h <- f.(as.matrix(Z[,1:existing.dishes]),i,lam=function(s,t) inv(s,t,D))
        h <- ifelse(Z[i,1:existing.dishes]==1,h,1-h)
      }
      p <- ifelse(log,sum(g,h,p),prod(g,h,p))
    }
  } else {
    p <- prod(dpois(0,a/(1:N)))
  }

  p
}

# Just a comment:
# Simulate Properties:
#Zs <- foreach(i=1:B) %dopar% raibp(N=2,a=a,D=D1)
B <- 1e5
a <- 1
Zs <- as.list(1:B); Zs <- lapply(Zs,function(x) {
                                      cat("\r: ",x/B)
                                      raibp(N=3,a=a,D=D1)
                                      })
EZ <- sum.matrices(Zs)/length(Zs)
EZ; apply(EZ,1,sum)
UZ <- unique.matrix(Zs)
L <- length(UZ[[2]])

# STOPPED HERE: WHAT's WRONG? TRY FOR r=1:10
get.one <- function(l,a,eps=1e-5) {
  B <- length(Zs)
  z <- UZ[[2]][[l]]
  p <- UZ[[1]][l] / B
  ci <- qnorm(c(.025,.975),p,sqrt(p*(1-p) / B))
  In <- F
  error <- F
  di <- dibp(z,a)
  di2<- dibp2(z,a)
  da <- daibp(z,a)
  if (abs(di-di2)>eps || abs(di-da)>eps) error <- T
  if (p > ci[1] && p < ci[2]) In <- T

  out <- c(p,dibp(z,a),dibp2(z,a),daibp(z,a),In,error)
  #names(out) <- c("Empirical","dibp","dibp2","daibp","IN","Error")
  out
}

result <- Rapply(as.list(1:L),function(x) get.one(x,a))
colnames(result) <- c("Empirical","dibp","dibp2","daibp","IN","Error")
apply(result,2,sum)

cumProb <- 0
for (i in 1:L) {
  cumProb <- cumProb + daibp(UZ[[2]][[i]],a)
}
cumProb

# aibp:
Zs <- as.list(1:B); Zs <- lapply(Zs,function(x) {
                                      cat("\r: ",x/B)
                                      raibp(N=3,a=a,D=D)
                                      })
UZ <- unique.matrix(Zs)
r <- 5;  Z <- UZ[[2]][[r]]; Z; 
UZ[[1]][r]/length(Zs); daibp(Z,a)

