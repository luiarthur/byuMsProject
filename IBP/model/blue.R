as.hex.string <- function(n){
hex.digits <- character(16);
for (i in 1:9) hex.digits[i] <- as.character(i);
for (i in 10:15) hex.digits[i] <- LETTERS[i-9];
hex.digits[16] <- as.character(0);
if (n<16) H <- paste(as.character(0),hex.digits[n],sep=character(1));
if (n==0) H <-paste(as.character(0), as.character(0), sep=character(1));
if (n>=16) {
H1 <- hex.digits[n%/%16];
n <- n-(n%/%16)*16;
H2 <- hex.digits[n];
if (n==0) H2 <- as.character(0);
H <- paste(H1,H2,sep=character(1));
}
H;
}

as.colour.string <- function(R,G,B){
r <- as.hex.string(R);
g <- as.hex.string(G);
b <- as.hex.string(B);
paste(rawToChar(as.raw(35)),r,g,b,sep=character(1));
}

X <- seq(0,255,1);
BLUE <- numeric(256);
#for (i in 1:256){BLUE[i] <- as.colour.string(0,0,X[i]);}
#M <- matrix(0,10,10);
#for (i in 1:10) for (j in 1:10) M[i,j] <- (sample(256,1)-1);
#image(M,col=BLUE,axes=FALSE)
