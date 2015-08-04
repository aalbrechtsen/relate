############ do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l)
  unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
Args<-function(l,args){
 if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
  cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
  q("no")
}
 arguments<-list()
 for(a in names(args))
   arguments[[a]]<-getArgs(a,l)

 if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
   cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
   q("no")
 }
 for(a in names(args))
   if(is.null(arguments[[a]]))
     arguments[[a]]<-args[[match(a,names(args))]]

   
 arguments
}

print.args<-function(args,des=as.list("")){
  cat("needed arguments:\n")
  mapply(function(x)cat("->  ",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("optional arguments (defaults):\n")
  mapply(function(x)cat("->  ",x," (",args[[x]],") :",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}
###### ####### ###### ###### ###### #######
# choose your parameters and defaults
# NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(file="linkage.pdf",stat=NULL,perm=NULL,chr=NA,pos=NA,alpha="0.05",ylab="Statistic",xlab="SNP position",main="relatedness mapping",w=10,h=6,type="l")

#optional description list to be passed to print.args - does not need to contain all arguments
des<-list(file="name of the pdf file",stat="the output filename with the statistic (stat.x)",perm="the output filename with the permutation results (x.txt)",chr="the filename of the chromosomes (stripped.chr)",pos="the filename of the positions (stripped.pos)",alpha="the significance threshold",ylab="label for the y-axis",xlab="label for the x-axis",main="label of plot",w="width of the pdf in cm",h="height of the pdf in cm",type="type of points (p,l,b)")

#if no argument aree given prints the need arguments and the optional ones with default
if(length(l)==0) print.args(args,des)
#get arguments and add to workspace
attach(Args(l,args))




#####################################
#### add your code here
stat<-scan(stat)
perm<-scan(perm)
alpha<-as.numeric(unlist(strsplit(alpha,",")))
sigThres<-quantile(perm,1-alpha,type=1)
cat("The significanse threshold is",sigThres,"with an alpha of",alpha,"\n")


pdf(file,w=as.numeric(w),h=as.numeric(h))
{
  if(is.na(pos)){
    pos<-1:length(stat)
    plot(pos,stat,xlab=xlab,ylab=ylab,main=main,type=type)

  }
  else{
    pos<-scan(pos)
    if(!is.na(chr)){
      chr<-scan(chr)
      m<-tapply(pos,chr,max)
      pos<- c(0,cumsum(m))[chr]+pos
    }
    p<-plot(pos,stat,xlab=xlab,ylab=ylab,main=main,type=type)
    if(exists("m")){
      abline(v=m[-length(m)])
      text(m/2+cumsum(c(0,m[-length(m)])),rep(mean(range(stat)),length(m)),names(m))
    } 
  
  }
  abline(h=sigThres,lty=2)
  text(rep(pos[round(length(pos)/40)]),sigThres-diff(range(stat))/30,alpha)
}
dev.off()


