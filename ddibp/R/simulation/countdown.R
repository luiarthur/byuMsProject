count.down <- function(old.time,i,B) {
  prog <- round(100*i/B,4)
  new.time <- Sys.time()
  time.diff <- as.numeric(new.time-old.time)
  time.remain <- time.diff * (B-i)
  if (time.remain < 60) {
    secs <- round(time.remain)
    time.remain <- paste0(secs,"s     ")
  } else if (time.remain<3600) {
    mins <- round(time.remain%/%60)
    secs <- round(time.remain%%60)
    time.remain <- paste0(mins,"m ",secs,"s        ")
  } else {
    hrs <- round(time.remain%/%3600)
    mins <- round((time.remain%%3600) %/% 60)
    time.remain <- paste0(hrs,"h ",mins,"m         ")
  }
  cat(paste0("\rProgress: ",prog,"%. Time Remaining: ",time.remain," "))
  if (i==B) cat("100%\n")
}

a.image <- function(Q,color=rev(heat.colors(100)),#paste0("gray",100:0),
                    numbers=F,num.cex=1,
                    numcolor="black",axis.num=T,...) {

  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
  
  ec <- apply(Q,2,sum) 
  er <- apply(Q,1,sum)
  seq1 <- seq(0,1,len=length(ec))
  seq2 <- seq(1,0,len=length(er))

  if (axis.num) {
    axis(1,at=seq1,lab=ec)
    axis(2,at=seq2,lab=er,las=2,...)
  }

  if (numbers) {
    xx <- rep(1:ncol(Q),each=nrow(Q))
    yy <- rep(1:nrow(Q),ncol(Q))
    text(seq1[xx],seq2[yy],c(Q),col=numcolor,font=2,cex=num.cex)
    #print(t(Q)[xx,yy])
    #for (x in 1:ncol(Q)) {
    #  for (y in 1:nrow(Q)) {
    #    text(seq1[x],seq2[y],t(Q)[x,y],col=numcolor,font=2,cex=num.cex)
    #  }
    #}  
  }
}
