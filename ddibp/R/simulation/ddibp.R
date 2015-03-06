source("countdown.R")
source("aibp2.R")
a.image <- function(Q,color=rev(heat.colors(100)),#paste0("gray",100:0),
                    numbers=F,num.cex=1,
                    numcolor="black",axis.num=T,...) {

  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
  
  ec <- apply(Q,2,sum) 
  er <- apply(Q,1,sum)
  seq1 <- seq(0,1,len=length(ec))
  seq2 <- seq(1,0,len=length(er))

  if (axis.num) {
    axis(1,at=seq1,lab=ec)
    axis(2,at=seq2,lab=er,las=2,...)
  }

  if (numbers) {
    xx <- rep(1:ncol(Q),each=nrow(Q))
    yy <- rep(1:nrow(Q),ncol(Q))
    text(seq1[xx],seq2[yy],c(Q),col=numcolor,font=2,cex=num.cex)
    #print(t(Q)[xx,yy])
    #for (x in 1:ncol(Q)) {
    #  for (y in 1:nrow(Q)) {
    #    text(seq1[x],seq2[y],t(Q)[x,y],col=numcolor,font=2,cex=num.cex)
    #  }
    #}  
  }
}

#X <- as.matrix(iris[1:5,1:4])
#X <- matrix(c(1,2,3,2,3,4,3,4,5),3,3,byrow=T)


rddibp <- function(n=3,a=1,D=NULL,fun=function(d) exp(-d)) {

  f.mat <- NULL
  if (is.null(D)) {
    f.mat <- matrix(1,n,n)
  } else {
    n <- nrow(D)
    f.mat <- fun(D)
  }
  f.mat[which(upper.tri(f.mat))] <- 0

  h <- apply(f.mat,1,sum)
  len <- length(h)
  S <- matrix(rep(h,rep(len,len)),len,len,byrow=T)
  A <- f.mat / S

  lam <- rpois(n,a/h) 
  #lam <- c(1,3,0,2,1)
  K <- sum(lam)
  set.K <- as.list(rep(0,n))
  if (lam[1]>0) set.K[[1]] <- 1:lam[1]
  mx <- 0
  if (n>1) {
    for (i in 2:n) {
      mx <- max(mx,max(set.K[[i-1]]))
      if (lam[i]>0) set.K[[i]] <- (mx+1):(mx+lam[i])
    }
  }

  get.owner.of.dish <- function(x) { # x is the dish number
    ind <- 1
    while (!(x %in% set.K[[ind]])) {
      ind <- ind + 1
    }
    ind # return the customer number of the owner of dish x
  }
  
  C <- matrix(0,n,K)
  # Is this right???
  if (K>0) {
    for (i in 1:n) {
      for (k in 1:K) {
        j <- sample(1:n,1,prob=A[i,])
        C[i,k] <- j
      }
    }
  }
  
  trace.inheritance <- function(cust,dish,visited=NULL) {
     #print(paste((cust,dish,i)))
     if (dish %in% set.K[[cust]]) {
       TRUE
     } else if (cust %in% visited) {
       FALSE
     } else {
       trace.inheritance(C[cust,dish],dish,c(visited,cust))
     }
  }

  #set.K; C
  #trace.inheritance(cust=1,dish=4)
  #trace.inheritance(cust=2,dish=4)
  #trace.inheritance(cust=3,dish=4)
  #trace.inheritance(cust=4,dish=4)
  #trace.inheritance(cust=5,dish=4)

  Z <- matrix(0,n,K)
  inherit.K <- as.list(1:n)
  if (K>0) {
    #for (k in 1:K) {
    #  for (i in 1:n) {
    #    Z[i,k] <- trace.inheritance(i,k)
    #  }
    #}
    for (i in 1:n) {
      Z[i,] <- apply(matrix(1:K),1,function(k) trace.inheritance(i,k))
    }
  }

  Z
}

