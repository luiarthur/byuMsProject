#rm(list=ls())
#options("width"=180)
#set.seed(1)
library(doMC); registerDoMC(strtoi(system("nproc",intern=T))/2)

# Plotting the Grid
#rIBP <- function(ppl=50,a=10,plotting=F){
#
#  last <- 0
#  #while(last==0) last <- rpois(1,a) #V1
#  last <- rpois(1,a)
#
#  #M <- matrix(0,ppl,last) #V1
#  M <- matrix(0,ppl,last) #V2
#  M[1,0:last] <- 1 #V1
#
#  for (n in 2:ppl){
#    for (k in 0:last){
#      mk <- sum(M[,k])
#      M[n,k] <- sample(0:1,1,prob=c(1-mk/n,mk/n)) #V1
#    }  
#    newLast <- last+rpois(1,a/n)
#    col0 <- matrix(0,ppl,newLast-last)
#    if (ncol(col0) > 0){ # added & newLast > 0
#      M <- cbind(M,col0)
#      M[n,(last+1):newLast] <- 1
#      last <- newLast
#    }
#  }
#
#  if (plotting) {
#    library(pheatmap)
#    pheatmap(M,cluster_row=F,cluster_col=F,display_numbers=T,number_format="%.0f")
#  }
#
#  M
#}


unique.matrix <- function(MM) {
  # Counts the number of unique matrices in a list.
  # MM = a list of matrices. We want the output to be:
  # 1) a list of UNIQUE matrices
  # 2) a vector of their counts

  eq.matrix <- function(A,B) {
    matrices.are.identical <- F
    if (all(dim(A) == dim(B))) {

      if((A == B) || dim(A)[2] == 0) matrices.are.identical <- T
    }
    matrices.are.identical
  }
  
  N <- length(MM)
  vec <- 1:N
  uniq.M <- list()
  
  counts <- c()
  taken <- rep(F,N)

  for (i in 1:N) {
    if (!taken[i]) {
      for (j in i:N) {
        if (!taken[j]) {
          if (eq.matrix(MM[[i]],MM[[j]])) {
            if (i==j) {
              counts <- c(counts,1)
              uniq.M[[length(uniq.M)+1]] <- MM[[j]]
            } else {
              counts[length(counts)] <- counts[length(counts)] + 1         
              taken[j] <- T
            }
          }
        }
      }
      taken[i] <- T
    }
  }
  list(counts,uniq.M)
}

#This was done on 20 September. I need to see if P[lof(Z)] = empirical probability.
my.K <- 3
draw.one <- function(ppl=3,a=1,MaxK=20){
  p <- rbeta(MaxK,a/MaxK,1)
  Z <- matrix(NA,nrow=ppl,ncol=MaxK)
  for ( k in 1:ncol(Z) ) {
    Z[,k] <- rbinom(ppl,1,p[k])
  }
  Z
}


#source("../model/lof.R",chdir=T)
#lof.Zs <- foreach(i=1:10000) %dopar% lof(rIBP(3,1))
##ZZ <- foreach(i=1:10000) %dopar% draw.one(3,1,20)
##lof.Zs <- lapply(ZZ)
#result <- unique.matrix(lof.Zs)
#result

# Also in ../model/model.R
#find.kh <- function(z){
#
#  # After finding kh,
#  # foreach unique vector,
#  # compute the factorial of the counts,
#  # then multiply the quantities together.
#  # So the constant in 14 is the product of the factorial of the counts of unique
#
#  #n <- nrow(z)
#  k <- ncol(z)
#  kh <- NULL
#  z.col   <- list()
#  for (i in 1:k) z.col[[i]] <- z[,i]
#  uniq.zk <- unique(z.col)
#  for (i in 1:length(uniq.zk)) kh[i] <- count(uniq.zk[[i]],z.col) 
#  kh
#}
#
#
## Also in ../model/model.R
#count <- function(m,Zs){ # counts the appearance of the matrix m in a list
#                         # of matrices Zs
#  cnt <- 0
#  for (i in 1:length(Zs)){
#    if (all(dim(Zs[[i]]) == dim(m))){
#      if (all(Zs[[i]] == m)) {
#        cnt <- cnt + 1
#      }
#    }
#  }
#  cnt
#}

#source("../model/model.R",chdir=T)
#plof <- function(Z,a=1,log=F) {
#  K <- ncol(Z)
#  n <- nrow(Z)
#  Kh <- find.kh(Z) # need this from ../model/model.R
#  Hn <- sum(1/(1:n))
#  mk <- apply(Z,2,sum)
#  
#  out <- NULL
#  if (!log) {
#    out <- a^K / prod(gamma(Kh+1)) * exp(-a*Hn) * prod(gamma(mk)*gamma(n-mk+1)/gamma(n+1))
#  } else {
#    out <- K * log(a) - sum(lgamma(Kh+1)) - a*Hn + sum(lgamma(mk)+lgamma(n-mk+1) - lgamma(n+1))
#  }
#
#  out
#}
#
#result[[1]]
#exp(plof(result[[2]][[9]],log=T))
#
#apply(matrix(1:length(result[[2]])),1,function(x) exp(plof(result[[2]][[x]],log=T)))


## Check to see if dIBP = empirial of P(Z)
#N <- 1e5
#n <- 2
#a <- .0001
#Zs <- foreach(i=1:N) %dopar% rIBP(n,a)
#result <- unique.matrix(Zs)
#result
#
#theory.prob <- apply(matrix(1:length(result[[2]])),1,function(x) dIBP(result[[2]][[x]],a,exch=F,const=T,log=F))
#
##cbind(normalize(theory.prob),normalize(result[[1]]))
#cbind(theory.prob,normalize(result[[1]]))
#na <- mean(unlist(lapply(Zs,sum))) #= Expected Number of Entries in Z, should be n*a
#
#apply(cbind(theory.prob,normalize(result[[1]])),2,sum) # Want: 1,1
#sum(unlist(lapply(result[[2]],sum)) * theory.prob)     # Want: n * a. Theory NOT ok.
#na                                                     # Want: n * a. Empirical ok.
#mean(unlist(lapply(Zs,ncol)))                          # Want: a*sum(1/(1:n)). OK

#source("../model/sumMatrices.R")
#summed.matrices <- sum.matrices(Zs,T)
#matrices <- sum.matrices(Zs,T)
#Some Assignment: ################################################
#
#alpha <- c(1,3,5,7)
#n <- 55
#M <- list(); length(M) <- length(alpha)
#counts <- list(); length(counts) <- length(alpha)
#props <- list(); length(props) <- length(alpha)
#results <- list(); length(results) <- length(alpha)

# Plots:
#png("buffet.png")
#  par(mfrow=c(2,2)) # Can't print plot?
#  for (i in 1:length(alpha)) {
#    M[[i]] <- plotGrid(n,alpha[i])
#
#    results[[i]] <- matrix(c(alpha[i],apply(M[[i]],2,sum),
#                           NA,apply(M[[i]],2,mean)), nrow=2, byrow=T)
#    colnames(results[[i]]) <- c("a",paste("k",1:ncol(M[[i]]),sep=""))
#    rownames(results[[i]]) <- c("Counts","Proportions")
#    results[[i]] <- t(results[[i]])
#  }
#  par(mfrow=c(1,1))
#dev.off()
#
## Format:
#for(i in 1:length(alpha)){
#  rows <- max(unlist(lapply(results,nrow))) - nrow(results[[i]])
#  zeros <- matrix(NA,rows,2)
#  results[[i]] <- rbind(results[[i]], zeros)
#}
#
#formatted.results <- cbind(results[[1]],results[[2]],results[[3]],results[[4]])
#rownames(formatted.results) <- c("a",paste("k",1:(nrow(formatted.results)-1),sep=''))
#formatted.results


