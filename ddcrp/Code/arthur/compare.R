# Progress:
pb <- function(i,n) { cat(paste0("\rProgress: ",round(i*1000/n)/10,"%")); if (i==n) cat("\n") } # arthur

source("../ddcrp-inference.R",chdir=T)
source("../data-modeling.R")
source("../../../code/auxGibbs/scala/crp/crp.R")

library(foreach) 
library(doMC)
registerDoMC(16)


alpha <- 3
N <- 10000
dat <- as.matrix(c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78))

do.it <- function(data=dat,a=alpha,n=N) {
  res <- opt.ddcrp.gibbs(dat=data, dist.fn=seq.dist, 
                     alpha=a, niter= n,
                     decay.fn=identity.decay,
                     lhood.fn= function(x) 0,
                     clust.traj=T,
                     summary.fn = ncomp.summary)

  make.uniq.lab <- function(vec=res$map.state[,"cluster"]) {
    uniqueLab <- unique(vec)
    for (i in 1:length(uniqueLab)) {
      ind <- vec==uniqueLab[i]
      vec[ind] <- i
    }
    vec
  }

  out <- NULL
  M <- res$clust

  if (!is.matrix(M)) {
    out <- matrix(make.uniq.lab(),1)
  } else {  
    out <- t(apply(M,1,make.uniq.lab))
  } 

  out
} 

compare <- function(niter=N,a=alpha) {
  blei.crp <- do.it(n=niter,a=a)
  blei.num.clust <- apply(blei.crp,1,max)
  dist.blei <- table(blei.num.clust) / sum(table(blei.num.clust))

  do.me <- function(i,B=niter) {pb(i,B); crp(n=length(dat),a=a)}
  my.crp <- foreach(i=1:niter,.combine=rbind) %dopar% do.me(i)
  my.num.clust <- apply(my.crp,1,max)
  dist.me <- table(my.num.clust) / sum(table(my.num.clust))
  
  list(alpha=a,blei.crp=blei.crp,my.crp=my.crp,
       dist.blei=dist.blei,dist.me=dist.me)
}

# Main: 
alphas <- as.list(c(1,3,5))
result <- lapply(alphas,function(x) compare(N,x))

sink("out/results.txt")
  O<-lapply(result, function(x){cat(paste0("alpha = ",x$a,"\n\n"))
                                print(x$dist.me);cat("\n")
                                print(x$dist.blei);cat("\n\n\n")})
sink()

