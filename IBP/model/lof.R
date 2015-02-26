lof <- function(Z) {
  z.col <- NULL
  #if (ncol(Z) <= 1) as.matrix(Z)
  Z <- as.matrix(Z)
  #if (ncol(Z) == 0) print(Z) 
  for (j in 0:ncol(Z)){
    z.col[j] <- paste(Z[,j],collapse="")
  }
  z.col.bin <- strtoi(z.col,base=2)
  
  lof.Z <- as.matrix(Z[,rev(order(z.col.bin))])
  while (ncol(lof.Z)>0 & sum(lof.Z[,ncol(lof.Z)])==0) 
     lof.Z <- as.matrix(lof.Z[,-ncol(lof.Z)])

  lof.Z
}

#library(pheatmap)
#Z <- matrix(sample(0:1,30,replace=T),5,6)
#lof.Z <- lof(Z)
#
#pdf("bayes/latex/report1/pics/Z.pdf")
#  pheatmap(Z,cluster_row=F,cluster_col=F,color=rev(gray.colors(100)),legend=F)
#dev.off()
#pdf("bayes/latex/report1/pics/lofZ.pdf")
#  pheatmap(lof.Z,cluster_row=F,cluster_col=F,color=rev(gray.colors(100)),legend=F)
#dev.off()
