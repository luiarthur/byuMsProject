# FUNCTIONS: ########################################################
Rapply <- function(L,f) { # L is a list, f is a function to apply to L[[x]]
                          # L apply takes a list, applies a function, 
                          # and rbinds it. Assumes output is vector.
  n <- length(L)
  out <- apply(matrix(1:n),1,function(i) f(L[[i]]))
  t(out)
}

unique.matrix <- function(X) {
  # Counts the number of unique matrices in a list.
  # X = a list of matrices. We want the output to be:
  # 1) a list of UNIQUE matrices
  # 2) a vector of their counts
  toMat <- function(s) {
    dims <- regexpr(": \\d* \\d*",s)
    begin <- as.integer(dims)+2
    end <- begin+attr(dims,"match.length")
    dims <- substr(s,begin,end)
    pos <- as.integer(regexpr(" ",dims))
    dims <- c(substr(dims,1,pos-1),substr(dims,pos+1,nchar(dims)))
    dims <- as.integer(dims)

    mat <- substr(s,1,begin-3)
    M <- matrix(0,dims[1],dims[2])
    if (mat>" ") {
      vec <- as.integer(strsplit(mat,",")[[1]])
      M <- matrix(vec,dims[1],dims[2])
    }

    M
  }

  S <- lapply(X,function(x) paste(toString(x),":",nrow(x),ncol(x),collapse=","))
  tab <- table(unlist(S))
  counts <- as.integer(tab)
  mat <- names(tab)
  uniq.M <- lapply(as.list(mat),toMat)
  
  ind <- sort(counts,index.return=T,decr=T)
  ind <- ind$ix
  list("count"=counts[ind],"matrix"=uniq.M[ind])
}

sum.matrices <- function(Ms,return.matrices=F) { 
# Ms is a list of matrices of different lengths
# return.matrices is a boolean. If FALSE, function returns the sum of the matrices.
# If TRUE, function returns a list of the matrices also.

  l <- length(Ms)
  max.c <- max(unlist(lapply(Ms,ncol)))
  max.r <- max(unlist(lapply(Ms,nrow)))
  
  for (i in 1:l) {
    M <- Ms[[i]]
    
    ncol0 <- max.c - ncol(M)
    nrow0 <- max.r - nrow(M)

    if (ncol0>0) {
      col0 <- matrix(0,nrow(M),ncol0)
      M <- Ms[[i]] <- cbind(Ms[[i]],col0)
    }

    if (nrow0>0) {
      row0 <- matrix(0,nrow0,ncol(M))
      M <- Ms[[i]] <- rbind(Ms[[i]],row0)
    }
  }

  out <- Reduce("+",Ms)

  if (return.matrices) out <- list("sum"=out,"matrices"=Ms)

  out
}
#####################################################################

get.new.dish <- function(z) {
  N <- nrow(z)
  K <- ncol(z)
  x <- rep(0,N)
  x[1] <- sum(z[1,])
  for (i in 2:N) {
    if (sum(x[1:(i-1)])+1 <= K) x[i] <- sum(z[i,(sum(x[1:(i-1)])+1):K])
  }
  x
}

inv <- function(s,t,d=D) 1/d[s,t] # inverse distance metric

# Calculates Probability of Customer_i Getting Dish_k
f. <- function(x,i=2,lam=inv,log=F) {
  K <- ncol(x)
  if (is.null(K)) K <- 0
  f <- function(x,i,k) {
    out <- 0
    if (sum(x[1:(i-1),k]>0)) {
      ind <- which(x[1:(i-1),k]==1)
      out <- sum(lam(ind,i))
    }
    out
  }

  out <- 0
  if (K>0) {
    t <- apply(matrix(1:K),1,function(k) f(x,i,k)) 
    out <- rep(0,K)
    if (sum(t)>0) {
      if (!log) {
        out <- t/sum(t)*sum(x[1:(i-1),])/i
        #out <- t^x[i,]*(sum(t)-t)^(i-x[i,])/sum(t) * sum(x[1:(i-1),])/i
      } else {
        #out <- x[i,]*log(t)+(i-x[i,])*log(sum(t)-t)-
        #       log(sum(t))+log(sum(x[1:(i-1),]))-log(i)
      }
    }
  }

  out
}


# For a GIVEN PERMUTATION!!!
raibp <- function(N=3,a=3,D=NULL,l=inv) {
  K <- rpois(1,a)
  Z <- matrix(0,N,K) 
  Z[1,0:K] <- 1 # The first customer draws a POI(a) number of new dishes
  
  # If no distance matrix is provided, customers will be equidistant.
  if (is.null(D)) {
    D <- matrix(1,N,N)
    #diag(D) <- 0
  }

  if (N>=2) {
    for (i in 2:N) {
      P <- f.(Z,i,lam=function(s,t,d) l(s,t,D))
      if (K>0) Z[i,] <- P > runif(K)

      newK <- K+rpois(1,a/i)
      col0 <- matrix(0,N,newK-K)

      if (ncol(col0) > 0) {
        Z <- cbind(Z,col0)
        Z[i,(K+1):newK] <- 1
        K <- newK
      }
    }
  }
  
  Z
}

ribp <- function(N=3,a=1){
  K <- rpois(1,a)
  Z <- matrix(0,N,K)
  Z[1,0:K] <- 1
  
  for (i in 2:N) {
    # Sample previously sampled dishes
    if (K>0) {
      mk <- apply(Z,2,sum)
      pk <- mk / i
      Z[i,] <- pk > runif(K)
    } 

    # Adding new dishes for customer i
    newK <- K+rpois(1,a/i)
    col0 <- matrix(0,N,newK-K)
    if (ncol(col0) > 0) {
      Z <- cbind(Z,col0)
      Z[i,(K+1):newK] <- 1
      K <- newK
    }
  }

  Z
}

dibp <- function(Z,a) {
  N <- nrow(Z)
  K <- ncol(Z)
  Hn <- sum(1/1:N)
  x <- get.new.dish(Z)
  mk <- apply(Z,2,sum)
  
  p <- a^K / prod(gamma(x+1)) * exp(-a*Hn) *
       prod(gamma(N-mk+1)*gamma(mk)/gamma(N+1))
  p
}

dibp2 <- function(Z,a) {
  N <- nrow(Z)
  K <- ncol(Z)
  Hn <- sum(1/1:N)
  x <- get.new.dish(Z)
  p <- a^K / prod(gamma(x+1)) * exp(-a*Hn)

  if (K>0) {
    for (i in 2:N) {
      p <- p / i^x[i]
      y <- sum(x[1:(i-1)])
      if (y>0) {
        for (k in 1:y) {# For each existing dish previous to customer i
          m_ik <- sum(Z[1:(i-1),k])
          p <- p * m_ik^Z[i,k] * (i-m_ik)^(1-Z[i,k]) / i
        }
      }
    }
  }

  p
}

daibp <- function(Z,a=3,D=NULL,l=inv,log=F) {
  N <- nrow(Z)
  K <- ncol(Z)
  x <- get.new.dish(Z)
  
  if (is.null(D)) {
    D <- matrix(1,N,N)
    diag(D) <- 0
  }

  p <- ifelse(log,0,1)
  if (K>0) {
    for (i in 1:N) {
      g <- dpois(x[i],a/i,log=log)
      h <- ifelse(log,0,1)
      existing.dishes <- sum(x[1:(i-1)])
      if (i>1 && existing.dishes>0) {
        h <- f.(as.matrix(Z[,1:existing.dishes]),i,lam=function(s,t) inv(s,t,D),log=log)
        h <- ifelse(Z[i,1:existing.dishes]==1,h,1-h)
      }
      p <- ifelse(log,sum(g,h,p),prod(g,h,p))
    }
  } else {
    p <- prod(dpois(0,a/(1:N)))
  }

  p
}
