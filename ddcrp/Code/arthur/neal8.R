# Progress:
pb <- function(i,n) { cat(paste0("\rProgress: ",round(i*1000/n)/10,"%")); if (i==n) cat("\n") } # arthur

library(foreach)
library(doMC)
registerDoMC(16)


y <- c(-1.48,-1.4,-1.16,-1.08,-1.02,.14,.51,.53,.78)
n <- length(y)


relabel <- function(x){
  uniq <- unique(x)
  y <- NULL
  for (i in 1:length(uniq)){
    y[x==uniq[i]] <- i
  }
  y
}

rg <- function(n=1){
  rnorm(n,0,1)
}

rF <- function(n=1,x){
  rnorm(n,x,.1)
}

dF <- function(x,phi){
  dnorm(x,phi,.1)
}

# Take out the c[i]
# c[i] needs to go back to one of the existing clusters
# OR it needs to get a new cluster.
# In practice, a=1, m=1.

c <- rep(1,n) 
phi <- rg(length(unique(c))+1)

# Chinese Restaurant:
crp <- function() {
  for (i in 1:n){ # start from two because the first cluster is 1
    k <- length(unique(c[-i]))
    h <- k + 1
    #phi <- phi[1:h]
    c[-i] <- relabel(c[-i])

    if (any(c[-i]==c[i])){
      phi[h] <- rg(1) # when m = 1
    }
    
    w <- NULL; w[h] <- 1
    for (t in 1:k){
      w[t] <- sum(c[-i]==t)
    }

    samp <- 1:h
    prob <- w * dF(y[i],phi) #a=1, m=1
    c[i] <- sample(samp,1,prob=prob)
  }
  c
}

do.neal <- function(i,B) { pb(i,B); neal <- crp() }

# Main:
B = 100000
neal.crp <- foreach(i=1:B,.combine=rbind,.errorhandling="remove") %dopar% do.neal(i,B)
neal.num.clust <- apply(neal.crp,1,max)
dist.neal <- table(neal.num.clust) / sum(table(neal.num.clust))
round(dist.neal,4)

