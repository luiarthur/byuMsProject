source("blue.R")
source("gibbs.ibp.R") # Gibbs sampler
source("model.R")
source("sumMatrices.R")
source("countdown.R")
library(corrplot)
#a.image <- function(Q,color=BLUE,...) {
a.image <- function(Q,color=paste0("gray",0:100),...) {
  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
}

a1 <- matrix(c(0,1,0,0,0,0,
               1,1,1,0,0,0,
               0,1,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0),6,6,byrow=T)

a2 <- matrix(c(0,0,0,1,1,1,
               0,0,0,1,0,1,
               0,0,0,1,1,1,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0),6,6,byrow=T)

a3 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               1,0,0,0,0,0,
               1,1,0,0,0,0,
               1,1,1,0,0,0),6,6,byrow=T)

a4 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,1,1,1,
               0,0,0,0,1,0,
               0,0,0,0,1,0),6,6,byrow=T)

# New Stuff. Generate Data.
N <- 100
K <- 4
Z <- matrix(sample(0:1,N*K,replace=T),N,K)
for (i in 1:N) {
  while (sum(Z[i,])==0) {
    Z[i,] <- sample(0:1,K,replace=T)
  }  
}
A <- matrix(c(a1,a2,a3,a4),K,byrow=T)
ZA <- Z%*%A
sigX <- .5
E <- matrix(rnorm(prod(dim(ZA)),0,sigX),dim(ZA)[1],dim(ZA)[2])
Y <- ZA + E


gibbs.post <- function(X=Y,siga=1,sigx=sigX,a=1,B=1000,burn=B*.1,showProgress=T,
                       a.a=1,a.b=1,plotProgress=F) {
  B <- ceiling(B/50)*50 

  D <- ncol(X)
  N <- nrow(X)

  Z <- rIBP(N,.5)
  alpha <- NULL
  alpha[1] <- a

  Hn <- sum(1/(1:N))

  p.x.z <- function(Z,log=T,with.norm.const=T) { # p.x.z = Likelihood
    K <- ncol(Z)
    H <- Z %*% solve(t(Z)%*%Z+(sigx/siga)^2 * diag(K)) %*% t(Z)
    I <- diag(nrow(H))

    if (!log) {
      norm.const <- ifelse(with.norm.const,
                           (2*pi)^(N*D/2) * sigx^((N-K)*D) * siga^(K*D),1)
      out <- 
      ( norm.const * det(t(Z)%*%Z + (sigx/siga)^2 * diag(K))^(D/2) )^(-1) *
      exp(
        -1/(2*sigx^2) * tr(t(X) %*% (I-H) %*% X)
      )
    } else {
      norm.const <- ifelse(with.norm.const,
                           (N*D/2)*log(2*pi)+((N-K)*D)*log(sigx)+(K*D)*log(siga),0)
      out <- 
      -( norm.const + (D/2)*log(det(t(Z)%*%Z + (sigx/siga)^2 * diag(K))) ) + 
      -1/(2*sigx^2) * tr(t(X) %*% (I-H) %*% X)
    }
    #out + dIBP(Z,a,log=T,exch=F,const=F)
    out
  } # p.x.z
  
  #p.x.z(z) # Remove this line!!!  
  p.zik.x <- function(z,i,k,log=T) { # P[z_{ik}=1|z_{-i,k}].  # exact

    zz <- z
    zz[i,k] <- 1
    a <- p.x.z(zz,log=T,with.norm.const=F) 
    zz[i,k] <- 0
    b <- p.x.z(zz,log=T,with.norm.const=F) 
    mk = sum(z[-i,k])
    p1 <- mk/N # p = Prior
    p0 <- 1 - p1

    p0 <- b+log(p0)#7Nov,2014  
    p1 <- a+log(p1)#7Nov,2014
    
    if (mk==0){
      out <- 0 # essentially doing nothing. because mk=0 => z[i,k] = 0.
    } else if (!log) { # Not Logged
      #out <- a * p / (a*p + b*(1-p))
      out <- 1 / (1+exp(p0-p1))
    } else { # Logged
      #print(paste("a:",a,"b:",b))
      #out <- a + log(p) - log(taylor.e(a)*p + taylor.e(b)*(1-p)) # HERE IS MY PROBLEM! 31 OCT
      #out <- a + log(p) - log(exp(a)*p + exp(b)*(1-p)) # HERE IS MY PROBLEM! 31 OCT
      out <- -log(1+exp(p0-p1))
    }
    
    out
  }
  
  # DOES THIS MAKE SENSE???
  # STOPPED HERE 8 NOVEMBER. START AGAIN on 10 NOVEMBER.

  # Need a way to sample alpha. Draw posterior alpha.
  #sampNewAlpha <- function(alpha=a,z,s=5) {
  #  prior <- function(x) dgamma(x,1,1,log=T)
  #  lp <- prior(seq(0,5,len=100))
  #  ll <- p.x.z(z)
  #  
  #}

  sampNewCols <- function(z,i,s=9,a=alpha[1]) { #s is the max num of columns. Avoid mh.
    prior <- function(l) dpois(l,a/N,log=T)  

    lp <- prior(0:s) # log prior
    ll <- apply(matrix(0:s),1,function(x) { 
                                col1 <- matrix(0,N,x)
                                col1[i,] <- 1
                                p.x.z(cbind(z,col1))
                                }) # log like
    lg <- lp+ll

    lpi <- apply(matrix(0:s),1,function(i) -log(sum(exp(lg-lg[i+1]))))
    Pi <- exp(lpi)
    
    #print(Pi)
    sample(0:s,1,prob=Pi)
  }

  Zs <- as.list(1:B)
  Zs[[1]] <- matrix(1,N)

  for (b in 2:B) { # B = num of iterations in Gibbs
    old.time <- Sys.time()
    z <- Zs[[b-1]]
    alpha[b] <- alpha[b-1]

    for (i in 1:N) { # iterate through all rows of Z
      K <- ncol(z)

      k <- 1
      while(k<=K) { # iterate through all columns of Z
        if (K>0) {
          #p <- p.zik.x(z,i,k,log=T) # pzik=1|x. Here
          p <- p.zik.x(z,i,k,log=F) # pzik=1|x.
          #cat(paste("\r Num of Z columns:",K))
          #Sys.sleep(.01)
          #u <- log(runif(1))
          u <- runif(1)
          if (p<=0) { # HERE!!!!!!!!!!!!!!!!!!!!!!!
            #if (p==0) print("Removed a column")
            z <- as.matrix(z[,-k])
            k <- k-1
          } else if (p>u){
            z[i,k] <- 1
          } else {
            z[i,k] <- 0
          }
        }

        k <- k + 1
        K <- ncol(z)
      } # end of its through all columns of Z

      # Drawing new dishes should be prior times likelihood
      #new <- rpois(1,a/N)
      new <- sampNewCols(z,i,a=alpha[b])

      if (new>0) {
        #print(paste("a/N = ", a/N ,". New Columns:",new,"ncol(z) =",ncol(z)))
        col1 <- matrix(0,N,new)
        col1[i,] <- 1
        z <- as.matrix(cbind(z,col1))
      }
    } # end of its through all rows of Z

    # Z|alpha propto alpha^{a-1} exp{-alpha Hn}
    #   alpha propto alpha^{a-1} exp{-alpha / b}
    # if a=1, b=1
    # a|Z ~ Gamma(a+K,(1/b+Hn)^(-1)) 
    #a.a <- a.a+ncol(z)
    #a.b <- (1/a.b+Hn)^(-1)
    #alpha[b] <- rgamma(1,a.a,scale=a.b)
    alpha[b] <- rgamma(1,a.a+ncol(z),scale=1/(1/a.b+Hn))
    Zs[[b]] <- z
    
    if (b %% 50 == 0) {
      gp <- b%/%50
      sink("draw.post.out/Z.post.results",append=T)
        Zs[ ((gp-1)*50+1) : (gp*50) ]
      sink()  
    } else if (b==1) {
      sink("draw.post.out/Z.post.results")
      sink()  
    }

    if (plotProgress) {
      n.col <- unlist(lapply(Zs[1:b],ncol))
      plot(n.col,xlab="Iteration",ylab="K+",
           main=paste0("Columns of Z ","(",b,")"),col="pink",lwd=3,type="b",pch=20)
      abline(h=mean(n.col),col="blue",lwd=3)     
    }

    if (showProgress) count.down(old.time,b,B)#pb(b,B)
  } # end of gibbs

  #out <- Zs[(burn+1):B]
  out <- list("Zs"=Zs,"alpha"=alpha)
  out
}

#R <- as.matrix(iris[,-5])
#f <- iris[,5]
#elapsed.time <- system.time(M <- gibbs.post(R,a=2,B=1000,burn=0,showProgress=T))

elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=1000,burn=0,showProgress=T,
                                              plotProgress=T))
M <- out$Zs
alpha <- out$alpha

sink("draw.post.out/Z.post.results")
  M
  elapsed.time
sink()  

burn <- 100#round(length(M) * .1)

sink("draw.post.out/Z.matrix")
  Z
sink()  

sink("draw.post.out/A.matrix")
  A
sink()

sink("draw.post.out/Y.matrix")
  Y
sink()

n.col <- unlist(lapply(M,ncol))
pdf("draw.post.out/traceplot.pdf")
  plot(n.col,type="l",main="Trace Plot: Number of Columns in Z",lwd=1,cex=.1,
       col="blue",pch=20)
  abline(h=mean(n.col[-(1:burn)]),lwd=2,col="red")
dev.off()


EAXZ <- function(X,Z,siga=1,sigx=sigX) {
  k <- ncol(Z)
  Ik <- diag(k)
  ZT <- t(Z)
  out <- solve(ZT%*%Z +(sigx/siga)^2 * Ik,ZT%*%X)
  #out <- solve(ZT%*%Z, ZT%*%X)
  out
}

Z.post <- M[-(1:burn)] # Burn in about 100
Z.post.mean <- sum.matrices(Z.post) / length(Z.post)
#Z.post.mean <- ifelse(Z.post.mean>runif(length(Z.post.mean)),1,0)
Z.post.mean <- ifelse(Z.post.mean>.9,1,0)
a.image(Z.post.mean)


col0.ind <- which(apply(Z.post.mean,2,function(x) sum(x)==0))
Z.post.mean <- Z.post.mean[,-col0.ind]
#Z.post.mean <- cbind(Z.post.mean[,3],Z.post.mean[,2],
#                     Z.post.mean[,4],Z.post.mean[,1])



one.A <- EAXZ(Y,Z.post.mean,siga=1,sigx=sigX)
d2 <- 2
d1 <- ceiling(nrow(one.A)/d2)


pdf("draw.post.out/postA.pdf")
  a.image(one.A,main="Posterior Mean for A")
dev.off()

plot.post.As <- function() {
  par(mfrow=c(d1,d2))
  for (i in 1:nrow(one.A)) {
    one.Ai <- matrix(one.A[i,],6,6) # matrix(Y[n,],6,6) = X[[n]]
    a.image(one.Ai,main=paste0("Posterior Mean A",i))
    #a.image(one.Ai,main=paste0("Posterior Mean A",i),col=BLUE)
  }
  par(mfrow=c(1,1))
}

pdf("draw.post.out/postA66.pdf")
  plot.post.As()
dev.off()

pdf("draw.post.out/Y.pdf")
  a.image(Y,main="Y")
dev.off()

pdf("draw.post.out/Z.pdf")
  a.image(Z,main="Z")
dev.off()

pdf("draw.post.out/A.pdf")
  a.image(A,main="A")
dev.off()

pdf("draw.post.out/postZ.pdf")
  a.image(Z.post.mean,main="Posterior Estimate for Z")
dev.off()

plot.As <- function() {
  par(mfrow=c(2,2))
    a.image(a1,main="A1")
    a.image(a2,main="A2")
    a.image(a3,main="A3")
    a.image(a4,main="A4")
  par(mfrow=c(1,1))
}

pdf("draw.post.out/a6by6.pdf")
  plot.As()
dev.off()


check.accuracy <- function(i) {
  sum(Z[i,]) == sum(Z.post.mean[i,])
}
num.acc <- apply(matrix(1:nrow(Z.post.mean)),1,check.accuracy)
mean(num.acc)
# mean(Z == Z.post.mean) = .98

post.ZA <- Z.post.mean %*% one.A
plot.post.ZA <- function(n) {
  par(mfrow=c(1,2))
    a.image(matrix(Y[n,],6,6),main=paste0("n=",n,":  ",toString(Z[n,])))
    a.image(matrix(post.ZA[n,],6,6),
            main=paste0("n=",n,":  ",toString(Z.post.mean[n,])))
  par(mfrow=c(1,1))
}

X11(); plot.post.As()
X11(); plot.As()
X11()
plot.post.ZA(31)

#z <- sum.matrices(Z.post) / length(Z.post)
#z.i <- matrix(z[6,],1)
#x.i <- matrix(X[[6]],1)
#classes[[6]]
#eaxz <- EAXZ(x.i,z.i)
#image(t(matrix(eaxz[4,],6,6)),col=gray.colors(36))

#Z.post.mean <- sum.matrices(Z.post) / length(Z.post)
#Z.post.mean[7,]
#classes[[7]]

#mean.one.A <- apply(one.A,2,mean)
#mean.one.A <- matrix(mean.one.A,6,6) # matrix(Y[n,],6,6) = X[[n]]
#corrplot(stand(mean.one.A),"square")

#library(ggplot2)
#image(t(Z.post[[n]]))
#head(Z.post[[n]])
#apply(Z.post[[800]],2,mean)

#sink("temp.txt") 
#  Z.post
#sink()

