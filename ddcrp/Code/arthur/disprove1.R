source("../ddcrp-inference.R",chdir=T)
source("../data-modeling.R")

sim <- function(x) {
  
  #D <- matrix(c(0,1,.2,.3,
  #              1,0,x,4,
  #              .2,x,0,.5,
  #              .3,4,.5,0),4,4,byrow=T)

  x <- runif(6,0,2)
  D <- matrix(c(0,x[1],x[2],x[3],
                x[1],0,x[4],x[5],
                x[2],x[4],0,x[6],
                x[3],x[5],x[6],0),4,4,byrow=T)
  

  (D <- R[[2]])
  blei.ddcrp <-  opt.ddcrp.gibbs(dat=matrix(1,4,1), dist.fn=matrix.dist.fn(D), 
                                 alpha=1, niter= 1000,
                                 decay.fn=exp.decay(1),
                                 lhood.fn= function(x) 0,
                                 clust.traj=T,
                                 summary.fn = ncomp.summary)

  blei.traj <- blei.ddcrp$clust.traj

  d23 <- D[1,2]
  p23 <- mean(blei.traj[,1]==blei.traj[,2])
  d24 <- D[1,3]
  p24 <- mean(blei.traj[,1]==blei.traj[,3])
  matrix(c(d23,d24,p23,p24),2,2)


  if ((d23<d24 & p23<p24) || (d23>d24 & p23>p24)) {
    D
  } else {
    NULL
  }  
}

library(doMC);registerDoMC(16)
res <- foreach(i=1:100) %dopar% sim(i)
R <- Filter(Negate(is.null),res)

