# Next steps: Work on setting a log parameter
source("aibp.R")

# Setup:
N <- 5
a <- 1
d.lower <- 1
d.upper <- 100
D <- matrix(0,N,N)
D[lower.tri(D)] <- D[upper.tri(D)] <- round(runif(N*(N-1)/2,d.lower,d.upper))
B <- 1e4

# Simulation:
# AIBP:
  Zs <- lapply(as.list(1:B),function(x) {cat("\r",x/B); raibp(N=N,a=a,D=D)})
  UZ <- unique.matrix(Zs)
  EZ <- sum.matrices(Zs) / B
  apply(EZ,1,sum)
  apply(EZ,2,sum)
  UZ$matrix[1:5]
  daibp(UZ$matrix[[1]],a,D)
  daibp(UZ$matrix[[400]],a,D)
# AIBP with D equidistant:
  zs <- lapply(as.list(1:B),function(x) {cat("\r",x/B); raibp(N=N,a=a)})
  uz <- unique.matrix(zs)
  ez <- sum.matrices(zs) / B
  apply(ez,1,sum)
  apply(ez,2,sum)
  uz$matrix[1:5]
  daibp(UZ$matrix[[1]],a)
  daibp(UZ$matrix[[400]],a)
# IBP:
  zi <- lapply(as.list(1:B),function(x) {cat("\r",x/B); ribp(N=N,a=a)})
  ui <- unique.matrix(zi)
  ei <- sum.matrices(zi) / B
  apply(ei,1,sum)
  apply(ei,2,sum)
  ui$matrix[1:5]
  dibp(UZ$matrix[[1]],a)
  dibp(UZ$matrix[[400]],a)
# Distribution of counts seems to behave similarly. 
# The corresponding matrices are different.
# Distribution of row sums is unchanged.
# IS Distribution of column sums is unchanged?
# But surely, given a previous configuration (before i), i will change. Check.
  UZ$count[1:100]
  uz$count[1:100]
  ui$count[1:100]
  UZ$matrix[[8]]
  uz$matrix[[8]]
  ui$matrix[[8]]

# Not Logged:
aibp.empiric <- UZ$count / B
aibp.theory  <- unlist(lapply(UZ$matrix,function(x) {daibp(x,a,D)}))
aibp.result  <- cbind(aibp.empiric,aibp.theory)
aibp.result[50:60,]
tail(aibp.result)
# Logged:
#laibp.empiric <- log(aibp.empiric)
#laibp.theory  <- unlist(lapply(UZ$matrix,function(x) daibp(x,a,D,log=T)))
#laibp.result  <- cbind(laibp.empiric,laibp.theory)
#head(laibp.result)

