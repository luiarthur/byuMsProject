source("../model.R",chdir=T)
library(pheatmap)

n <- 5
d <- 6
m <- matrix(0,n,d)

B <- 10000
registerDoMC(10)
M.list <- foreach(i=1:B) %dopar% rmatnorm(m,(1:nrow(m))*diag(nrow(m)),diag(ncol(m)))
EM <- Reduce("+",M.list) / B

pheatmap(EM,cluster_row=F,cluster_col=F,
         display_numbers=T,number_format="%.2f",
         color=gray.colors(100,start=.7,end=.3),
         main="Matrix Normal(0,1,1)")

