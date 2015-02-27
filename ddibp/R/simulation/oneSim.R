B <- 1e4 # problems
a <- .5

# 1 close to 2, 2 close-ish to 3, 1 far from 3. 
D1 <- matrix(c(0,1,9,
               1,0,3,
               9,3,0),3,3)

# 1 far from to 2, 2 close-ish to 3, 1 close from 3. 
D2 <- matrix(c(0,9,1,
               9,0,3,
               1,3,0),3,3)

# 1 close to 2, 2 far from 3, 1 close-ish 3. 
D3 <- matrix(c(0,1,3,
               1,0,9,
               3,9,0),3,3)

#D <- matrix(0,3,3); D[which(lower.tri(D))] <- 1; 
#D <- matrix(1,3,3); D[which(upper.tri(D))] <- 0; 
#D <- matrix(c(0,1,20,1000,
#              1,0,50,30,
#              20,50,0,4,
#              1000,30,4,0),4,4)

exp.decay <- function(s,t,d) ifelse(s>t,0,exp(-d[s,t]))
inv <- function(s,t,d) ifelse(s>t,0,1/d[s,t])

one.sim <- function(D.name,a,B=1e4,num.cex=1) {
  D <- eval(parse(text=D.name))

  cat("Getting Draws (1/4): \n")
  Z.o <- lapply(as.list(1:B),function(x) {ot <- Sys.time() 
                                          o <- raibp(N=nrow(D),a=a)
                                          count.down(ot,x,B); o})
  Z.d <- lapply(as.list(1:B),function(x) {ot <- Sys.time(); o <- rddibp(a=a,D=D)
                                          count.down(ot,x,B); o}) 
  Z.a <- lapply(as.list(1:B),function(x) {ot <- Sys.time()
                                          o <- raibp(N=nrow(D),a=a,D=D,l=exp.decay)
                                          count.down(ot,x,B); o})
  #Z.a <- lapply(as.list(1:B),function(x) raibp(N=nrow(D),a=a,D=D,l=inv.decay))
  #Z.a <- lapply(as.list(1:B),function(x) raibp(N=nrow(D),a=a))

  cat("Calculating Expected Values (2/4): \n")
  EZO <- sum.matrices(Z.o)/B
  EZD <- sum.matrices(Z.d)/B
  EZA <- sum.matrices(Z.a)/B

  # Expected number of columns are different for ddIBP.
  mncolo <- mean(unlist(lapply(Z.o,ncol)))
  mncold <- mean(unlist(lapply(Z.d,ncol)))
  mncola <- mean(unlist(lapply(Z.a,ncol)))


  cat("Getting Unique Matrices from Draws (3/4): \n")
  uzo <- unique.matrix(Z.o)
  uzd <- unique.matrix(Z.d)
  uza <- unique.matrix(Z.a)

  # make a table:
  compare <- function() {

    options("width"=180)
    u.all <- unique.matrix(c(Z.o,Z.d,Z.a))
    n <- length(u.all$c)
    info <- matrix(0,n,7)
    colnames(info) <- c("matrix#",paste(rep(c("IBP","aIBP","ddIBP"),each=2),
                        c("empirical","theoretical")))

    cmo <- get.freqs(Z.o)
    cma <- get.freqs(Z.a)
    cmd <- get.freqs(Z.d)

    info[,2] <- unlist(lapply(u.all$m,function(z) get.freq(z,cmo)))
    info[,3] <- unlist(lapply(u.all$m,function(z) daibp(z,a=a)))
    info[,4] <- unlist(lapply(u.all$m,function(z) get.freq(z,cma)))
    info[,5] <- unlist(lapply(u.all$m,function(z) daibp(z,a=a,D=D,l=exp.decay)))
    info[,6] <- unlist(lapply(u.all$m,function(z) get.freq(z,cmd)))
    info[,7] <- NA
    info <- info[order(info[,3],decreasing=TRUE),]
    info[,1] <- 1:n

    list("info"=info,"unique"=u.all)
  }

  cat("Comparing Methods (4/4): \n")
  M <- compare()

  options("width"=80)
  list("EZO"=EZO,"EZD"=EZD,"EZA"=EZA,"uzo"=uzo,"uzd"=uzd,"uza"=uza,"M"=M,
       "mncolo"=mncolo,"mncola"=mncola,"mncold"=mncold)
}

