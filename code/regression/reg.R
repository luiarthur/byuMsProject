#The TWO functions below are inverses of each other:
  #Takes a DISTANCE MATRIX and returns the lower/upper triangle
  matrix.to.column <- function(M){
    M[which(lower.tri(M))]
  }

  #Takes a column and transforms back to the DISTANCE MATRIX
  column.to.matrix <- function(V){
    l <- length(V)
    h <- (-1 + sqrt(1+8*l)) / 2 + 1
    y <- matrix(0,h,h)
    y[lower.tri(y)] <- V
    y <- y + t(y)
    y
  }

#Takes a list of matrices and returns the weights for each matrix
get.columns <- function(Ms){
  c <- length(Ms)
  d <- dim(Ms[[1]])[1]
  r <- (1 + d) * d / 2
  M <- matrix(0,r,c)
  
  for (i in i:n){
    M[,i] <- matrix.to.column(Ms[[i]])
  }

  mod <- lm(M[,1] ~ M[,-1])

  return(mod)
  #mod$Coefficients returns the weights for each matrix
}

