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

unique.matrix <- function(X) {
  # Counts the number of unique matrices in a list.
  # X = a list of matrices. We want the output to be:
  # 1) a list of UNIQUE matrices
  # 2) a vector of their counts

  S <- lapply(X,function(x) paste(toString(x),":",nrow(x),ncol(x),collapse=","))
  tab <- table(unlist(S))
  counts <- as.integer(tab)
  mat <- names(tab)
  uniq.M <- lapply(as.list(mat),toMat)
  
  ind <- sort(counts,index.return=T,decr=T)
  ind <- ind$ix
  list(counts[ind],uniq.M[ind])
}

Rapply <- function(L,f) { # L is a list, f is a function to apply to L[[x]]
                          # L apply takes a list, applies a function, 
                          # and rbinds it. Assumes output is vector.
  n <- length(L)
  out <- apply(matrix(1:n),1,function(i) f(L[[i]]))
  t(out)
}

get.freq <- function(m,Zs) {
  N <- length(Zs)
  m.name <- paste(toString(x),":",nrow(m),ncol(m),collapse=",")

  S <- lapply(X,function(x) paste(toString(x),":",nrow(x),ncol(x),collapse=","))
  tab <- table(unlist(S))
  counts <- as.integer(tab)
  mat <- names(tab)

  count <- 0
  if (m.name %in% mat) {
    count <- counts[which(mat==m.name)]
  }

  count/N
}
