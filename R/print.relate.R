
print.relate<-function(x,...){
 if (!inherits(x, "relate")) 
        stop("Not an object of class HMMrelate!")
print(structure(list("k0,k1,k2"=round(rev(x$kResult),2),k.like=round(x$kLike,3),k.r=signif(x$kr,3),a=signif(x$a,3),u.like=round(x$uLike,1),po.like=round(x$poLike,1)),class="power.htest"),...)
}
