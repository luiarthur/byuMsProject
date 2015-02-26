source("../ddcrp-inference.R",chdir=T)
source("../data-modeling.R")

#D <- matrix(c(.0,.6,.2,
#              .6,.0,.1,
#              .2,.1,.0),3,3,byrow=T)


sim <- function(){
  
  x <- runif(3,0,.5)
  
  D <- matrix(c(0   ,x[1],x[2],
                x[1],0   ,x[3],
                x[2],x[3],0    ),3,3,byrow=T)

  (D <- R[[21]]) # 19
  blei.ddcrp <-  opt.ddcrp.gibbs(dat=matrix(1,3,1), dist.fn=matrix.dist.fn(D), 
                                 alpha=1, niter= 1000,
                                 decay.fn=exp.decay(1),
                                 lhood.fn= function(x) 0,
                                 clust.traj=T,
                                 summary.fn = ncomp.summary)

  blei.traj <- blei.ddcrp$clust.traj

  d12 <- D[1,2]
  d13 <- D[1,3]
  p12 <- mean(blei.traj[,1]==blei.traj[,2])
  p13 <- mean(blei.traj[,1]==blei.traj[,3])
  matrix(c(d12,p12,d13,p13),2,2,byrow=T);D 

  if ( (d12 < d13 && p12 < p13) || (d12 > d13 && p12 > p13) ) {
    D
  } else {
    NULL
  }
}

library(doMC);registerDoMC(16)
res <- foreach(i=1:100) %dopar% sim()
R <- Filter(Negate(is.null), res)

