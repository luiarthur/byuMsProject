#source("shallot/R/shallot.R")

library(jvmr)
library(shallot)
data <- c (a=-1.48,b=-1.4,c=-1.16,d=-1.08,e=-1.02,f=0.14,g=0.51,h=0.53,i=0.78)
lambda <- .1^-2
mu0 <- 0
lambda0 <- 1^-2

jainneal.log.predictive.density <- function(i,subset) {
  precision <- lambda0 + length(subset) * lambda
  mean <- (lambda0*mu0 +lambda*sum(data[subset])) / precision
  dnorm(data[i],mean,1/sqrt(lambda*precision/(precision+lambda)),log=T)
}

mcmc.params <- mcmc.parameters(log.density.NAME="jainneal.log.predictive.density")

dist <- ewens(mass(1,fixed=T),length(data))
mcmc <- collect(dist,n.draws=1000,mcmc.parameters=mcmc.params)
out <- process(mcmc)

