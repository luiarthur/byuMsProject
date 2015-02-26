tr <- function(M){
  sum(diag(M))
}


prob.X.given.Z.A.sigX <- function(X,Z,A,sigX){
  N <- nrow(X)
  D <- ncol(X)
  C <- X-Z%*%A
  (2*pi*sigX^2)^(-N*D/2) * exp(-1/(2*sigX^2) * tr(t(C)%*%C))
}


prob.A <- function(A,sigA){
  (2*pi*sigA^2)^(-K*D/2) * exp(-1/(2*sigA^2) * tr(t(A)%*%A))
}


prob.X.given.Z.sigX.sigA <- function(X,Z,sigX,sigA){
  N <- nrow(X)
  D <- ncol(X)
  K <- ncol(Z)
  Q <- t(Z) %*% Z + (sigX/sigA)^2 %*% diag(K)

  1 / ( (2*pi)^(N*D/2) * sigX^((N-K)*D) * sigA^(K*D) * det(Q)^(D/2) ) *
  exp(-1/(2*sigX^2) * tr(t(X) %*% (diag(N)-Z%*%solve(Q) %*% t(Z)) %*% X))
}
