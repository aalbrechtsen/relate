 
runHmmld.plink <- function(bed,bim,fam,pair=c(1,2),par=NULL,min.maf=0,LD="rsq2",epsilon=0.01,back=5,alim=c(0.001,0.15),start=NULL,prune=NULL,ld_adj=TRUE,fix.a=NULL,fix.k2=NULL,calc.a=FALSE,phi=0.013,timesToRun=10,timesToConverge=5,giveCrap=FALSE,convTol=0.1,method=1,back2=2*back) {
  #prepre plink stuff
  lb <- nchar(bed);
  ext <- substr(bed, lb-3, lb);
  if (ext==".bed") {
    stub <- substr(bed, 1, lb-4)
  } else {
    stub <- bed
    bed <- paste(bed, ".bed", sep="")
  }
  if (missing(bim))
    bim <- paste(stub, ".bim", sep="")
  if (missing(fam))
    fam <- paste(stub, ".fam", sep="")
  df.bim <- read.table(bim,as.is=T)
  snps <- as.character(df.bim[,2])
  if(giveCrap)
    cat("\t-> plinkTalking: Number of SNPs: ", length(snps), "\n")
  if (any(duplicated(snps)))
    stop("Duplicated SNP name(s)")
  df.fam <- read.table(fam)
  ped <- as.character(df.fam[,1])
  id <- as.character(df.fam[,2])
  if(giveCrap)
    cat("\t->plinkTalking: Number of individuals ",length(id),"\n length of snps",length(snps) )
  if (any(duplicated(id))) {
    cat("\t->PlinkTalking: Subject IDs not unique - concatenating with pedigree IDs\n")
    id <- paste(ped, id, sep=":")
    if (any(duplicated(id)))
      stop("Sample IDs are still not unique")
  }
  
  chr <- as.vector(t(df.bim[1]))
  position <- as.numeric(as.vector(t(df.bim[4])))/1000000
  keep <- as.integer(t(chr))
  keep[keep>22]<-0
  keep[is.na(keep)] <- 0
  keep<-as.logical(keep)
  cat("\t->plinkTalking: dataset contains ", sum(keep), "autosomale SNP's out of ",length(keep), "\n" )

  if(length(pair)!=2)
    stop("pair must contain the number of two individuals")

  if(pair[1]==pair[2])
    stop("pair must be two different individuals")

#  if(max(pair)>nrow(data)|min(pair)<1)
#    stop("pair must be integers between 1 and the number of rows in data")

  if(min.maf<0|min.maf>=1)
    stop("min.maf must have a value between 0 and 1")

  if(!LD%in%c("rsq2","D","D'"))
    stop("LD must be either \"rsq2\",\"D\" og \"D'\"")

  if(epsilon<0|epsilon>=1)
    stop("epsilon must have a value between 0 and 1")

  if(!is.null(prune)&back==0)
    stop("back must be greater than 0 when pruning")

  if(all(c(!is.null(fix.k2),fix.k2==0,calc.a))){
    cat("\tfixk2=0 and calc.a=T is a one-dimensional optimization will set\n\t\t timesToRun=3,timesToConverge=2")
    timesToRun=3
    timesToConverge=2
  }
  
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

  bedfullname <- file.path(dirname(bed),basename(bed))
  res<-.Call("doPlink",bedfullname,length(id),length(snps),as.integer(pair),as.numeric(position[keep]),pa,as.numeric(min.maf),as.integer(double_recom),as.integer(LD_var),as.numeric(epsilon),as.numeric(back),as.numeric(alim),as.integer(optim),as.numeric(start),prune,as.integer(ld_adj),fixA,fixK,as.integer(chr[keep]),as.integer(moment),calc.a,myphi,as.integer(timesToRun),as.integer(timesToConverge),as.integer(giveCrap),as.numeric(convTol),as.integer(method),as.integer(back2),as.integer(keep),PACKAGE="Relate")
  if (is.null(chr))
    res$chr=NULL
names(res$kResult)<-c("IBD=2","IBD=1","IBD=0")

  class(res)<-"HMMrelate"

  #res <- list(res,plinkKeep=keep,chromosomes=chr,positions=position)
  return(res)
}

