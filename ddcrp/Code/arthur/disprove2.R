source("../ddcrp-inference.R",chdir=T)
source("../data-modeling.R")

D <- matrix(c(.0,.1,.2,
              .1,.0,.3,
              .2,.3,.0),3,3,byrow=T)


blei.ddcrp <-  opt.ddcrp.gibbs(dat=matrix(1,3,1), dist.fn=matrix.dist.fn(D), 
                               alpha=1, niter= 1000,
                               decay.fn=exp.decay(1),
                               lhood.fn= function(x) 0,
                               clust.traj=T,
                               summary.fn = ncomp.summary)

blei.traj <- blei.ddcrp$clust.traj
mean(blei.traj[,2]==blei.traj[,1])
mean(blei.traj[,2]==blei.traj[,3])

