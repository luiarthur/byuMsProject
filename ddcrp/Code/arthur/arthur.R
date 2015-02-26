# Progress:
#pb <- function(i,n) { cat(paste0("\rProgress: ",round(i*1000/n)/10,"%")); if (i==n) cat("\n") } # arthur

source("../ddcrp-inference.R")
source("../data-modeling.R")


#docs <- read.documents("~/data/science/data/sci90/sci90-mult.dat")
#voc <- readLines("~/data/science/data/sci90/sci90-vocab.dat")
#dat <- corpus.to.matrix(docs[1:500], voc)
#
#res <- ddcrp.gibbs(dat=dat[1:100,], dist.fn=seq.dist, alpha=1,
#                   decay.fn=window.decay(100),
#                   doc.lhood.fn(0.5), 5, summary.fn = ncomp.summary)



dat <- as.matrix(c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78))
#dat <- as.matrix(sort(rnorm(100)))
plot(dat)
#n <- 30  
#dat1 <- sort(rnorm(n,0))
#dat2 <- sort(rnorm(n,10))
#dat3 <- sort(rnorm(n,20))
#dat <- as.matrix(c(dat1,dat2,dat3))
##dat <- (c(dat1,dat2,dat3))
#plot(1:n,dat1,xlim=c(0,3*n),ylim=c(min(dat),max(dat)),col="red",cex=3)
#points((n+1):(2*n),dat2,col="blue",cex=3)
#points((2*n+1):(3*n),dat3,col="green",cex=3)


res <- ddcrp.gibbs(dat=dat, dist.fn=seq.dist, # dist.fn=matrix.dist.fn(dist.matrix). See decay.R
                   alpha=1, niter= 10,
                   decay.fn=identity.decay,
                   #lhood.fn=function(x) exch.dirichlet.lhood(x,1),#function(x) sum(dnorm(x,10,2,log=T)),#? doc.lhood.fn(1), 
                   lhood.fn= function(x) 0,#doc.lhood.fn(1), 
                   summary.fn = ncomp.summary)
#res

#plot(dat)
uniqueLab <- unique(res$map.state[,"cluster"])
for (i in 1:length(uniqueLab)) {
  points(which(res$map.state[,"cluster"]==uniqueLab[i]),
     dat[which(res$map.state[,"cluster"]==uniqueLab[i])],
     col=rainbow(length(uniqueLab))[i],pch=20)
}
legend("topleft",legend=c(paste("Number of Clusters =",length(uniqueLab)),
                          paste("Number of Data Points =",length(dat))))

res
table(res$map.state[,"cluster"])

# Next Task:
# Find the Distribution of the Prior (i.e. since lhood= 0)

# Do Many Times:

do.it <- function() {
  res <- ddcrp.gibbs(dat=dat, dist.fn=seq.dist, # dist.fn=matrix.dist.fn(dist.matrix). See decay.R
                     alpha=1, niter= 10,
                     decay.fn=identity.decay,
                     #lhood.fn=function(x) exch.dirichlet.lhood(x,1),#function(x) sum(dnorm(x,10,2,log=T)),#? doc.lhood.fn(1), 
                     lhood.fn= function(x) 0,#doc.lhood.fn(1), 
                     summary.fn = ncomp.summary)
  table(res$map.state[,"cluster"])
}

library(foreach)
library(doMC)
registerDoMC(16)

N <- 1000
result <- foreach(1:N) %dopar% do.it()
len <- lapply(result,function(x) length(x))
EX <- Reduce("+",len) / N
EX
