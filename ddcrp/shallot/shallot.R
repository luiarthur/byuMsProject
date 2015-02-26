data <- c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
sigma <- 0.1
mu0 <- 0.0
sigma0 <- 1.0

s2 <- sigma*sigma
s02 <- sigma0*sigma0
s02Inv <- 1/s02


shallot.n <- length(data)

shallot.sample <- function(subset) {
  if ( all(is.na(subset)) ) rnorm(1,mean=mu0,sd=sigma0)
  else {
    s <- sum(data[subset])
    variance <- 1 / (s02Inv + length(subset) / s2)
    mean <- variance * (mu0 / s02 + s / s2)
    rnorm(1,mean=mean,sd=sqrt(variance))
  }
}

shallot.log.density <- function(i,parameter,subset) { # What is the subset param?
  dnorm(data[i],mean=parameter,sd=sigma,log=TRUE)
}


