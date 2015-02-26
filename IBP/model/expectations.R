#source("../IBP/IBP.R",chdir=T)
#options("width"=180)

# Plotting the Grid
rIBP <- function(ppl=50,a=10,plotting=F){

  last <- 0
  #while(last==0) last <- rpois(1,a) #V1
  last <- rpois(1,a) # I assume that the first customer may sample NO dishes

  #M <- matrix(0,ppl,last) #V1
  M <- matrix(0,ppl,last) #V2
  M[1,0:last] <- 1 #V1

  for (n in 2:ppl){
    for (k in 0:last){
      mk <- sum(M[,k])
      M[n,k] <- sample(0:1,1,prob=c(1-mk/n,mk/n))
    }  
    newLast <- last+rpois(1,a/n)
    col0 <- matrix(0,ppl,newLast-last)
    if (ncol(col0) > 0){ # added & newLast > 0
      M <- cbind(M,col0)
      M[n,(last+1):newLast] <- 1
      last <- newLast
    }
  }

  if (plotting) {
    library(pheatmap)
    pheatmap(M,cluster_row=F,cluster_col=F,display_numbers=T,number_format="%.0f")
  }

  M
}

expected.IBP <- function(n=20,a=7,N=1000){
  Z.list <- foreach(k=1:N) %dopar% rIBP(n,a,plotting=F) # Takes the longest time
  k.max <- function(z.list) max(unlist(lapply(z.list,ncol)))
  kmax <- k.max(Z.list)

  sum.each.list <- function(x.list,r=1,c=1){
    lower.bound <- function(one.matrix) ncol(one.matrix) >= c
    get.rc <- function(one.matrix) one.matrix[r,c]
    list.index <- which(unlist(lapply(x.list,lower.bound)))

    out <- sum(unlist(lapply(x.list[list.index],get.rc)))
    out
  }

  each.count <- matrix(0,n,kmax)
  for (r in 1:n){
    each.count[r,] <- foreach(c=1:kmax,.combine=cbind) %dopar% sum.each.list(Z.list,r,c)
  }

  rownames(each.count) <- paste("n",1:n,"=",apply(each.count,1,sum),sep="")
  colnames(each.count) <- paste("k",1:kmax,"=",apply(each.count,2,sum),sep="")

  each.count
}


library(pheatmap)
plot.expected.IBP <- function(M,a,disNum=F,prop=T,N,row.names=T,col.names=T,legend=T,
                              main=paste("Expected Draw from an IBP with a =",a)){
  nf <- "%.0f"
  if (prop) {
    nf <- "%.2f"
    erow <- apply(M,1,function(m) sum(m) / N)
    n <- nrow(M); k <- ncol(M)
    ecol <- apply(M,2,function(m) sum(m) / N / n)

    rownames(M) <- paste("n",1:n,"=",round(erow,5),sep="")
    colnames(M) <- paste("k",1:k,"=",round(ecol,5),sep="")
  }

  if (!row.names) rownames(M) <- NULL
  if (!col.names) colnames(M) <- NULL

  pheatmap(M,cluster_row=F,cluster_col=F,
           display_numbers=disNum,number_format=nf,
           #color=gray.colors(100,start=.7,end=.3),
           main=main,legend=legend)
}


EIBP <- function(a,makepdf=F,N=10000,prop=T) {
  M <- expected.IBP(n=20,a=a,N=N)
  filename <- paste("alpha",a,".pdf",sep="")
  if (makepdf){
    pdf(filename); plot.expected.IBP(M,a=a,prop=prop,N=N); dev.off()
  } else {
    plot.expected.IBP(M,a=a,prop=prop,N=N)
  }

  M
}

#alpha <- c(1,3,5,7)
#temp <- EIBP(1,F)
#plot.expected.IBP(temp,7,N=10000)
#result <- foreach(i=1:length(alpha)) %dopar% EIBP(alpha[i],T)
#a <- 7
#pdf("bayes/latex/report1/pics/draw1.pdf")
#  plot.expected.IBP(rIBP(10,a),a,N=1,main=paste("One Draw from IBP with a =",a),le=F)
#dev.off()
#pdf("bayes/latex/report1/pics/draw2.pdf")
#  plot.expected.IBP(rIBP(10,a),a,N=1,main=paste("One Draw from IBP with a =",a),le=F)
#dev.off()
#pdf("bayes/latex/report1/exp.pdf")
#  EIBP(a)
#dev.off()
################### Summary of Interesting Results: #########################
#
# 1) The proportion of customers that try dish k increases with k,
#                                is not affected for each k as a increases,
#                                decreases for each k as n increases.
#
# 2) The expected number of dishes taken by EACH customer ~= a.
# 3) The proportion of customers that take each dish is not clear. But seems
#    to follow a pattern.
#
#############################################################################

############################## Next Task: ###################################
#
# Changing Alpha:
#   1) E[X|Z=z]
#   2) E[X] = E[E[X|Z]]
#
# Metropolis:
#   1) X|Z ~ N(ZJ, s2_uI, s2I), where X(nxd), Z(nxk), J(kxd), I(dxd), I(kxk)
#   2) Z ~ IBP(a)
#
#   Q) What does X look like for various a, s2?
#   Q) What is Z|X?
#
#############################################################################

