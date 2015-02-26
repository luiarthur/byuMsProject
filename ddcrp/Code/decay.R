
# decay functions and distance functions

window.decay <- function(w)
  function (x) (as.integer(x <= w))

logistic.decay <- function(a, b=1)
  function (x) logistic((-x+a)/b)

exp.decay <- function(a, b = 1.0)
  function (x) exp(-x/a) * b + 1e-6

identity.decay <- function (x) # modified by Arthur Lui. 14 May, 2014.
{
  if (length(x) == 1) {
    if ( is.infinite(x) ) {
      0
    } else {
      1
    }  
  } else {
    #stop("I wasn't expecting this.")
    ifelse(is.infinite(x),0,1)
  }
}

euclidean.dist <- function(p1, p2)
{
  v <- p1 - p2
  sqrt(v %*% v)
}

manhattan.dist <- function(p1, p2)
{
  v <- abs((p1[1] - p2[1])) + abs((p1[2] - p2[2]))
}

subtract.dist <- function(p1, p2)
{
  #p1 - p1 # blei's code. is this supposed to be p1 - p2?
  p1 - p2 # Arthur Lui. 14 May, 2014.
}

seq.dist <- function(i,j)
{
  if (j <= i)
    i - j
  else
    Inf
}

input.based.dist.fn <- function(input, dist.fn)
{
  function (i,j) { dist.fn(input[i], input[j]) }
}

matrix.dist.fn <- function(dist.matrix)
{
  function (i,j) dist.matrix[i,j]
}

link.dist.fn <- function(adj)
{
  function (i, j) if (adj[i,j]==0) { Inf } else { 1 }
}
