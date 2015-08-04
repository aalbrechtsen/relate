 
relate<-function(data,position,pair=c(1,2),par=NULL,min.maf=0,LD="rsq2",epsilon=0.01,back=5,alim=c(0.001,0.15),start=NULL,prune=NULL,ld_adj=TRUE,fix.a=NULL,fix.k2=NULL,chr=NULL,calc.a=TRUE,phi=0.013,timesToRun=10,timesToConverge=5,giveCrap=FALSE,convTol=0.1,method=1,back2=2*back){




#error messages
  if(class(data)!="matrix")
    stop("data must be a matrix")

  snp<-dim(data)[2]
  if(snp!=ncol(data))
    stop("the number of positions does not match the number of SNPs (coloums) in data")

  if(length(pair)!=2)
    stop("pair must contain the number of two individuals")

  if(pair[1]==pair[2])
    stop("pair must be two different individuals")

  if(max(pair)>nrow(data)|min(pair)<1)
    stop("pair must be integers between 1 and the number of rows in data")

  if(min.maf<0|min.maf>=1)
    stop("min.maf must have a value between 0 and 1")

  if(!LD%in%c("rsq2","D","D'"))
    stop("LD must be either \"rsq2\",\"D\" og \"D'\"")

  if(epsilon<0|epsilon>=1)
    stop("epsilon must have a value between 0 and 1")

  if(!is.null(prune)&back==0)
    stop("back must be greater than 0 when pruning")

  if(all(c(!is.null(fix.k2),fix.k2==0,calc.a))){
    print("fixk2=0 and calc.a=T is a one-dimensional optimization will set timesToRun=3,timesToConverge=2")
    timesToRun=3
    timesToConverge=2
  }
    
  d<-as.integer(data)
  d[is.na(d)]<-0
  d<-as.integer(d)
  data<-matrix(d,ncol=snp)

  if(is.null(colnames(data)))
    colnames(data)<-paste("snp",1:dim(data)[2],sep="")
  
  #check cromo
  chromo <- c()
  if(is.null(chr))
    chromo <- rep(1,ncol(data))
  else
    chromo <- as.integer(chr)
  #check fixA
  fixA <- c()
  if(is.null(fix.a))
    fixA <- NULL
  else
    fixA <- as.numeric(fix.a)
  
  #check fixK
   fixK <- c()
  if(is.null(fix.k2))
    fixK <- NULL
  else
    fixK <- as.numeric(fix.k2)

  #check prune
  if(!is.null(prune)){
    if(prune==0)
      prune<-NULL
    else
      prune <- as.numeric(prune)
  }

  #point likelihood
  pa <- c()
  if(is.null(par))
    pa <- NULL
  else
    pa <- as.numeric(par)

  myphi <- c()
  if(is.null(phi))
    myphi <- NULL
  else
    myphi <- as.numeric(phi)

 
  #these are not used by the user
  #maf
  moment=FALSE
  double_recom=FALSE
  maf=NULL
  
  optim=1
  LD_var <- ifelse(LD=="rsq2",TRUE,FALSE)
  res<-.Call(interface,data,as.integer(pair),as.numeric(position),pa,as.numeric(min.maf),as.integer(double_recom),as.integer(LD_var),as.numeric(epsilon),as.numeric(back),as.numeric(alim),as.integer(optim),as.numeric(start),prune,as.integer(ld_adj),fixA,fixK,as.integer(chromo),as.integer(moment),calc.a,myphi,as.integer(timesToRun),as.integer(timesToConverge),as.integer(giveCrap),as.numeric(convTol),as.integer(method),as.integer(back2),PACKAGE="Relate")
  if (is.null(chr))
    res$chr=NULL
names(res$kResult)<-c("IBD=2","IBD=1","IBD=0")

  class(res)<-"relate"
  return(res)
}
