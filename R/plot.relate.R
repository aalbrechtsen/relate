 

plot.relate <- function(x,col=1:3,lwd=2,ylab="probability",xlab="position (Mb)",chr=NULL,legend=TRUE,...){
  post<-t(x$post)
  pos<-x$position
  if(is.null(x$chr)){
    plot(pos[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,ylab=ylab,xlab=xlab,...)
    lines(pos[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
    lines(pos[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
  }
  else{
    if(!is.null(chr)){
      pos<-pos[x$chr==chr]
      post<-post[x$chr==chr,]
      plot(pos[-length(pos)],post[-1,1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,ylab=ylab,xlab=xlab,...)
      lines(pos[-length(pos)],post[-1,2],col=col[2],lwd=lwd)
      lines(pos[-length(pos)],post[-1,3],col=col[3],lwd=lwd)
      
    }
    else{
      C<-names(table(x$chr))
      m<-c(0,cumsum(tapply(pos,x$chr,max)))
      pos2<-rep(NA,length(pos))
      for(tal in 1:length(C))
        pos2[x$chr==C[tal]]<-pos[x$chr==C[tal]]+m[tal]
      plot(pos2[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,ylab=ylab,xlab=xlab,...)
      lines(pos2[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
      lines(pos2[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
      abline(v=m)
      for(tal in 1:length(C))
        text(m[tal]+diff(m)[tal]/2,0.5,C[tal],col="gray")
    }
  }
  if(legend)
      legend(150,0.9,paste("IBD=",2:0),col=1:3,lty=1)

}
