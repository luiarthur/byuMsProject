rm(list=ls())
source("../model.R",chdir=T)
options("width"=60)
options("scipen"=10,digits=5)
#set.seed(1)

#(Z <- rIBP(ppl=3,a=1,plot=F)) # (Seeded) Generate a Z from IBP

#N <- 5
#X.list <- foreach(i=1:N) %dopar% 
#  r.X.given.Z(Z=Z,d=2,su=.1,sv=.1) # Generate X's given that Z

# Goal:
# - Recover the latent structure (Z) responsible
#   for generating the observed data, X.

log.lik <- function(x,z,a=a,su=1) {
  n <- nrow(x)
  k <- ncol(z)
  k1 <- find.k1(z)
  d <- ncol(x)
  j <- matrix(1,k,d)
  zj <- z%*%j
  Hn <- sum(1/(1:n))
  mk <- apply(z,2,sum)
  #A <- sum(log(factorial(n-mk)) + log(factorial(mk-1)))
  A <- sum(lgamma(n-mk+1) + lgamma(mk))

  #result <- -.5/su^2 * tr(t(x-zj)%*%(x-zj)) + k * log(a) -
  #          sum(log(factorial(k1))) - a*Hn - k*log(factorial(n)) - A
  result <- -.5/su^2 * tr(t(x-zj)%*%(x-zj)) 
  result
}

Gibbs <- function(x, init=rIBP(3,1,F), cs=list(diag(nrow(init)),diag(ncol(init))), 
                  a=1, N=10^4){

  out <- list(); length(out) <- N
  out[[1]] <- init
  cnt <- 0
  d <- ncol(x)
  n <- nrow(x)

  cat(paste(rep("#",50),collapse="")); cat("\n")
  for(i in 2:N){
    #print(i)
    z <- out[[i]] <- out[[i-1]]

    # Think of a proposal mechanism!!!
    #cand <- rmatnorm(z%*%J,cs[[1]],cs[[2]])
    #cand <- rmatnorm(z,cs[[1]],cs[[2]])
    # TRY THIS:
    #   Use small sample length(X.list) = 10'
    #   Proposoal = Prior => cand <- rIBP(a)
    #                        r <- ratio of likelihoods
    cand <- rIBP(n,a,F)

    r <- log.lik(x,cand,a) - log.lik(x,z,a)
    #if (!is.na(r) && r>log(runif(1))){
    if (r>log(runif(1))){
      out[[i]] <- cand
      cnt <- cnt + 1
    }

    if (i %% (N/50) == 0) cat(">")
  }
  cat("\n")
  print(cnt/N)
  return(out[-c(1:round(N/10))])
}

####################################### Tested my Prior: it's good
#N.x <- 5
#a <- 1
#Z <- rIBP(ppl=3,a=a,plot=F)
#X.list <- foreach(i=1:N.x) %dopar% 
#  r.X.given.Z(Z=Z,d=2,su=.1,sv=.1) # Generate X's given that Z
#N <- 10000
#Z.est <- Gibbs(X.list,a=a,N=N)
#
#uniq.Z.est <- unique(Z.est)
#est.uniq.counts <- as.vector(foreach(m=uniq.Z.est,.combine=rbind) %dopar% count(m,Z.est))
#est.uniq.props  <- est.uniq.counts / N 
#uniq.Z.est[which.max(est.uniq.props)]
#
#likely.ind <- which(est.uniq.props > .05)
#likely.Z <- uniq.Z.est[likely.ind]
#likely.p <- est.uniq.props[likely.ind]
#likely.Z

###########################################
#est.plot(uniq.props[uniq.props > .01],type='h')
#est.likely.Z <- uniq.Z.est[which(uniq.props > .01)]
#est.likely.Z.prop <- uniq.props[which(uniq.props > .01)]

#unique(Z.est)
#Z

# Test my Prior
##N <- 1000
##Zs <- foreach(i=1:N) %dopar% rIBP(3,.001,F)
##uniq.Z <- unique(Zs)
##
###counts <- as.vector(foreach(m=Zs,.combine=rbind) %dopar% count(m,Zs))
###props  <- counts / N 
##
##uniq.counts <- as.vector(foreach(m=uniq.Z,.combine=rbind) %dopar% count(m,Zs))
##uniq.props  <- uniq.counts / N 
##plot(uniq.props[uniq.props > .01],type='h')
##likely.Z <- uniq.Z[which(uniq.props > .01)]
##likely.Z.prop <- uniq.props[which(uniq.props > .01)]
##
##
##likely.Z.prop[12]
##dIBP(likely.Z[[12]],const=T)

#############################################################3
#N <- 100000
#a <- .1
#Zs <- foreach(i=1:N) %dopar% rIBP(3,a,F)
#uniq.Z <- unique(Zs)
#uniq.counts <- as.vector(foreach(m=uniq.Z,.combine=rbind) %dopar% count(m,Zs))
#uniq.props  <- uniq.counts / N  # The empirical prob of each unique Z
#f <- function(x) dIBP(x,a=a,exch=T) # T: Exchangeable => LOF
#uniq.dibp <- sapply(uniq.Z,f)  # The theory prob of each unique Z
#
#sum(uniq.props) # should be 1
#sum(uniq.dibp)  # should be 1
#
#uniq.props[1:10]
#uniq.dibp[1:10]
###norm.const <- (normalize(sapply(Zs,f)) / sapply(Zs,f))[1]
