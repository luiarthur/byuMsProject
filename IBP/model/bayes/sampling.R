# Goal: Retrieve Z and A given X
rm(list=ls())
source("../model.R",chdir=T)
source("../lof.R",chdir=T)
options("width"=200)

#set.seed(311) # March 11, 2014.
#set.seed(312) # March 11, 2014.

#Set Global Parameters for generating data###################################
n <- 100
a <- 2# a=1,3 - Good Results
kk <- 31
d <- 5

#############################################################################

genData <- function(Z,A,E=NULL) {
  # Function for generating observations X
  #   - Requires Feature Matrix Z, and Value Matrix A
  #   - Optional Error Matrix E
  if (is.null(E)) {
    n <- nrow(Z)
    k <- ncol(A)
    U <- diag(n)
    V <- diag(k)
    M <- matrix(0,n,k)
    E <- rmatnorm(M,U,V)
  }

  X <- Z%*%A+E
  X
}
  
llik <- function(x,z,A) {
  term <- x-z%*%A
  (-.5 * tr(t(term)%*%term) - .5*tr(t(A)%*%A))
}

llik.Z <- function(x,z,A) {
  term <- x-z%*%A
  -.5 * tr(t(term)%*%term) + dIBP(z,a,log=T)
}

llik.A <- function(x,z,A) {
  term <- x-z%*%A
  -.5 * tr(t(term)%*%term) - .5*tr(t(A)%*%A)/2
}

rqz <- function(z,span=.3) {
  z <- as.matrix(z)
  n <- nrow(z)
  k <- ncol(z)
  
  r <- runif(1) 

  changeDigits <- function(z) {
    samp <- sample(0:(n*k),round(n*k*span))
    z[samp] <- ifelse(z[samp] == 1,0,1)
    lof(z)
  }

  addCol <- function(z) {
    newZ <- cbind(z,sample(0:1,n,replace=T))
    lof(newZ)
  }

  rmCol <- function(z) {
    s <- sample(0:k,1) 
    newZ <- z[,-s]
    lof(newZ)
  }
  
  if (k>1){
    if(r < 1/3) {
      changeDigits(z)
    }  else if (r < 2/3) {
      addCol(z)
    } else {
      rmCol(z)
    }  
  } else if (k==1) {
    if (r < 1/2) {
      changeDigits(z)
    } else {
      addCol(z)
    }  
  } else {
    addCol(z)
  }

}


bayes <- function(x,a,k=kk,B=100000,spn=.5,...){ # a determines the number of features
  # Goal: 1) Choose (generate) Z and A
  #       2) Generate X_{nxd} = Z_{nxk}A_{kxd} + e_{nxd}, e ~ MN(0,I_n,I_d)
  #       3) Estimate (recover) latent Z and A given X

 
  # Initialize Acceptance Rate Counters
  cnt.Z <- cnt.A <- 0

  add0 <- function(z){
    # Function that adds appropriate number of columns of 0's to z because
    # z changes through each iteration. And for conformality, z needs to have
    # a fixed length.
    z <- cbind(z,matrix(0,nrow(z),k-ncol(z)))
    z
  }


  # Initialize output list:
  out <- list()   
  length(out) <- B
  ################


  # Picking a Z and A to generate data. 
  z <- add0(rIBP(n,a))
  A <- round(rmatnorm(matrix(5,k,d),1*diag(k),1*diag(d)))

 
  # Initialize first list element of output
  out[[1]] <- list("Z"=z,"A"=A)

  # Define d and n for easier reading of code
  d <- ncol(x)
  n <- nrow(x)
  
  # Metropolis
  # Progress bar. Just for display. Not important to computation
  #cat(paste(rep("#",50),collapse="")); cat("\n")
  pb <- txtProgressBar(2,B,width=30,...)

  for (i in 2:B){ # B is a big number

    D <- out[[i]] <- out[[i-1]]

    #Update Z:#######################################################################
    cand.Z <- add0(rqz(D$Z,span=spn)) # This isn't doing so well. about 60% to 70% accuracy.
    #cand.Z <- add0(rIBP(n,a,F)) # Using the prior as the proposal. 
               # Higher accuracy. Even when N = 100. about 80%.

    #r <- llik(x,cand.Z,D$A) - llik(x,D$Z,D$A) # An alternative way of implementation
    r <- llik.Z(x,cand.Z,D$A) - llik.Z(x,D$Z,D$A) # + logP(x->x') - logP(x'->x)
    if (r > log(runif(1))) {
      out[[i]]$Z <- cand.Z
      cnt.Z <- cnt.Z + 1
    }
    D <- out[[i]]
    #End of Update Z ################################################################


    #Update A:#######################################################################
    cand.A <- rmatnorm(D$A,1*diag(k),1*diag(d)) # Symmetric Proposal
    #cand.A <- rmatnorm(matrix(sample(5,k*d,replace=T),k,d),2*diag(k),1*diag(d))
    #r <- llik(x,t$Z,cand.A) - llik(x,t$Z,t$A)
    r <- llik.A(x,D$Z,cand.A) - llik.A(x,D$Z,D$A)
    if (r > log(runif(1))) {
      out[[i]]$A <- cand.A
      cnt.A <- cnt.A + 1
    }
    D <- out[[i]]
    #End of Update A ################################################################

    setTxtProgressBar(pb, i)
    #if(i%%(B/50)==0) cat(">") # Just another part of the progress bar
  }

  #cat("\n")
  close(pb)
  print(paste("Z Acceptance Rate: ",cnt.Z/B))
  print(paste("A Acceptance Rate: ",cnt.A/B))
  #Sys.sleep(2)
  #out[-c(1:(B%/%10))]
  out
}

#MAIN: ###################################################
one.sim <- function(it="",bayes.loop=100000,span=.5){ # 4 minutes for one run of bayes(B = 100000)
  Z <- rIBP(n,a) 
  while (ncol(Z) == 0) Z <- rIBP(n,a) # Just want Z matrices with >= 1 cols
  #while (ncol(Z) != 3) Z <- rIBP(n,a) # Just want Z matrices with 3 columns
  k <- ncol(Z)
  A <- round(rmatnorm(matrix(5,k,d),2*diag(k),2*diag(d)))
  X <- genData(Z,A)

  out <- bayes(X,a,B=bayes.loop,style=3,spn=span)
  B <- length(out)

  out.Z <- lapply(out,function(x) x$Z)
  mean.Z <- Reduce("+",out.Z) / B   # Version 2
  roundZ <- ifelse(mean.Z >=.5,1,0) # Version 2
  mean.Z <- as.matrix(mean.Z[,which(apply(mean.Z,2,sum)>0)]) # Trunc Column of 0's
  roundZ <- as.matrix(roundZ[,which(apply(roundZ,2,sum)>0)]) # Trunc Column of 0's
  add0 <- function(z) cbind(z,matrix(0,nrow(z),max(ncol(roundZ),ncol(Z))-ncol(z)))
  z.acc <- sum(add0(Z) == add0(roundZ)) / prod(dim(add0(Z)))
  

  out.A <- lapply(out,function(x) x$A)
  mean.A <- Reduce("+",out.A) / B
  length.uniq.A <- length(unique(out.A))
  cut.mean.A <- mean.A[0:ncol(roundZ),]
  #mean.A

  post.ZA <- roundZ%*%cut.mean.A
  
  result <- list("z.acc"=z.acc,"Z"=Z,"round.Z"=roundZ,"mean.Z"=mean.Z,"A"=A,
                 "mean.A"=round(cut.mean.A),"X"=X,"post.ZA"=post.ZA,"out"=out)

  sink(paste("results/result",it,".txt",sep="")); print(result); sink()  
  result
}

#temp <- one.sim("temp",1000)
registerDoMC(16)
N <- 16
#system.time(result <- foreach(i=1:N,.errorhandling="pass") %dopar% one.sim(i,100000))
#system.time(result <- foreach(i=1:N,.errorhandling="pass") %dopar% one.sim(i,1000))

# Diagnostics:
  out <- one.sim(bay=1000,span=.05)
  outs <- out$out
  Zs <- lapply(outs,function(x) x$Z)
  
  sums <- sapply(Zs,function(z) sum(z)); plot(sums,type='l',lwd=3)
