n <- k <- a <- d <- su <- sv <- Z <- J <- 0

if (set==1) {
# 1: This one is interesting ####################
n <- 5; k <- 2; a <- 5; d <- 5; su <- 7; sv <- 2#
Z <- matrix(c(1,0,0,1,1,1,0,1,1,0),n,byrow=T)   #  
J <- matrix(1,k,d)                              #
#################################################
} else if (set==2) {
# 2: This one is interesting ####################
n <- 5; k <- 2; a <- 5; d <- 5; su <- 7; sv <- 2#
Z <- matrix(c(0,0,0,0,1,1,0,0,0,0),n,byrow=T)   #
J <- matrix(1,k,d)                              #
#################################################
} else if (set==3) {
# 2: This one is interesting ####################
n <- 5; k <- 5; a <- 5; d <- 5; su <- 7; sv <- 2#
Z <- matrix(c(rep(c(1,1,1,.5,.5),5)),n,byrow=T)   #
J <- matrix(1,k,d)                              #
#################################################
} else if (set==4) {
# 2: This one is interesting ####################
n <- 5; k <- 5; a <- 5.5; d <- 5; su <- 7; sv <- 2#
Z <- matrix(c(rep(c(1,1,1,0,0),5)),n,byrow=T)   #
J <- matrix(1,k,d)                              #
#################################################
} else if (set==5) {
# 2: This one is interesting ####################
n <- 5; k <- 5; a <- 5.5; d <- 5; su <- 7; sv <- 2#
Z <- matrix(sample(0:1,n*k,replace=T),n,byrow=T)   #
J <- matrix(1,k,d)                              #
#################################################
} else if (set==6) {
# 2: This one is interesting ######################
n <- 3; a <- 5; d <- 2; su <- 7; sv <- 2          #
Z <- matrix(c(1,0,1,0,1,1),n)                     #
k <- ncol(Z)                                      #
J <- matrix(1,k,d)                                #
###################################################
} else if (set==7) {
# 2: This one is interesting ######################
n <- 3; a <- 5; d <- 2; su <- 7; sv <- 2          #
Z <- matrix(c(1,1,1,0,1,1),n)                     #
k <- ncol(Z)                                      #
J <- matrix(1,k,d)                                #
###################################################
} else if (set==8) {
# 2: This one is interesting ######################
n <- 3; a <- 5; d <- 2; su <- 7; sv <- 2          #
Z <- matrix(c(1,0,1,1,1,1),n)                     #
k <- ncol(Z)                                      #
J <- matrix(1,k,d)                                #
###################################################
} else if (set==9) {
# 2: This one is interesting ######################
n <- 3; a <- 5; d <- 2; su <- 7; sv <- 2          #
Z <- matrix(c(1,0,1,0,0,1),n)                     #
k <- ncol(Z)                                      #
J <- matrix(1,k,d)                                #
###################################################
} else if (set==10) {
# 2: This one is interesting ######################
n <- 3; a <- 5; d <- 2; su <- 7; sv <- 2          #
Z <- matrix(c(1,0,1,0,1,0),n)                     #
k <- ncol(Z)                                      #
J <- matrix(1,k,d)                                #
###################################################
} else if (set==11) {
# 2: This one is interesting ######################
n <- 3; a <- 5; d <- 2; su <- 7; sv <- 2          #
Z <- matrix(c(1,1,1,1,1,1),n)                     #
k <- ncol(Z)                                      #
J <- matrix(1,k,d)                                #
###################################################
}
