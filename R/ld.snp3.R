ld.snp3 <- function(snpdata, back=100){
  nams <- colnames(snpdata)
  snp<-dim(snpdata)[2]
  d<-as.integer(snpdata)
  d[is.na(d)]<-0
  d<-as.integer(d)
  snpdata<-matrix(d,ncol=snp)
  
  if(is.null(nams))
    colnames(snpdata)<-paste("snp",1:dim(snpdata)[2],sep="")
  else
    colnames(snpdata) <- nams
  ans <- .Call("snp_pair_ld", snpdata,back, PACKAGE = "Relate")
  ans
}


