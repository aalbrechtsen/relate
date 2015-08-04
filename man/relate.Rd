\name{relate}
\alias{relate}
\title{A function for estimating relatedness between to individuals using SNP data}
\description{
  This function estimates the probability of sharing alleles identity by descent (IBD) accross the genome. The method is based on a continouos time Markov model with hidden states. The hidden states are the IBD states between two diploid chromosome pairs. We assume that the individuals are not inbreed and thus the individuals can share 0, 1 or 2 alleles IBD. The SNPs are allowed to be in linkage disequilibrium (LD). To accomidate LD the methods need SNP for several individuas in order to estimate the allele frequencies and the pairwise LD. The method return the posterior probabilies of the IBD states accross the genome  and the overall IBD sharing. 
}
\usage{
relate(data,position,pair=c(1,2),par=NULL,min.maf=0,LD="rsq2",epsilon=0.01,back=5,alim=c(0.001,0.15),start=NULL,prune=NULL,ld_adj=TRUE,fix.a=NULL,fix.k2=NULL,chr=NULL,calc.a=TRUE,phi=0.013,timesToRun=10,timesToConverge=5,giveCrap=FALSE,convTol=0.1,method=1,back2=2*back)
}

\arguments{
  \item{data}{Integer Matrix. A  matrix with SNP genotypes where NA or 0 denotes missing data, 1 for AA, 2 for Aa  and 3 for aa. The number of individuals is the number of rows and the number of SNPs is the number of coloums}
  \item{position}{ The position of each SNP in centi Morgan (or mega bases). If centi Morgan is used then phi should be set to 1}
  \item{pair}{Integer vector of length two with the row numbers of the two individuals where relatedness is to be estimated}
  \item{par}{Optional numeric vector c(a,k2,k1,k0) of parameters used instead of optimazation}
  \item{min.maf}{The minumum minor allele frequecy allowed}
  \item{LD}{The measure used to select the previous SNP to condition on, ("D'","D" or "rsq2")}
  \item{epsilon}{The error rate}
  \item{back}{The number of privous SNPs that can be conditioned on (see details for recomadations) }
  \item{alim}{The allowed range for a}
  \item{start}{Optional starting point for the optimazition}
  \item{prune}{The maximum value allowed for pairwise LD. If 0 then to pruning is performed}
  \item{ld_adj}{Logical. use the pairwise emission probabilities to correct for LD}
  \item{fix.a}{Numeric. Fix the a parameters to this value}
  \item{fix.k2}{Numeric. Fix the k2 parameter to this value}
  \item{chr}{a vector of chromosomes numbers (only relevant when
    multiple chromosomes are used)}
  \item{calc.a}{Estimate the a parameter from the overall IBD sharing. appropreate for distantly related individuals or individual who are related though one or two paths of the same length}
  \item{phi}{Numeric. The recombination rate in Morgans per Mega base (m/Mb)}
  \item{timesToRun}{Integer. The maximum number of times the optimization is run}
  \item{timesToConverge}{Integer. The number of times the optimazation should reach the same optimum}
  \item{giveCrap}{Integer. If non zero then emission probabilies (given
  the unobserved state and allele phase) and the haplotype probabilities
  are returned. Also more runtime information is given.}
  \item{convTol}{Numeric. The tolerance for stating that the likelihood have reached the same likelihood.}
	\item{method}{Two version has been implemented, method=0 with back2=2*back should be the same as method=1. Method=0 is much faster when running allpairs}
	\item{back2}{This is mainly used for debuging}
}
\details{

How to select the number of privious SNP that can be conditioned on. First of all if is no LD in the data chosse back=0.If you want to prune SNP away based on LD and/or you want to accomidate LD in the model choose <back> as the number of SNP where you expect there to be LD. For example if you have 500,000 SNP you expect their to be LD between a lot of SNPs in a region. Here I use back=50. If you only have 10,000 SNPs I would use back=5. 


If there is LD in the data but you want to remove the LD before the analysis then set adj\_ld to FALSE, prune to some numeric value larger than zero (e.g. 0.2) and back to some number.
}
\value{
  \item{kResult}{The value for sharing 2, 1 or 0 alleles IBD}
  \item{kLike}{The maximal -log likelihood}
  \item{kr}{The co-ancestry coefficient}
  \item{a}{The rate of change between IBD states}
  \item{uLike}{The likehood for being unrelated}
  \item{LD}{The LD measure}
  \item{t}{distance between the used SNPs (some SNPs might have been discarded)}
  \item{snp}{The number of used SNPs (some SNPs might have been discarded)}
  \item{position}{The position of the used SNPs  (some SNPs might have been discarded)}
  \item{double_recom}{If false then no instanstationus double recombination is allowed}
  \item{alim}{the allowed range of a}
  \item{poLike}{The likelihood for being parent-offspring}
  \item{back}{The number of privious SNPs that have been conditioned on (and/or used for pruning)}
  \item{post}{a matrix with posterior probabilities for the hidden IBD states}
  \item{timesRun}{The number of times the optimazation algorihm have been run}
  \item{timesConverged}{The number of times the optimazation algorithm have found the same maximum}
  \item{convergenceInfo}{matrix with estimates of the overall IBD sharing and the likelihood for each of the optimaziations}
  \item{usedSnps}{a vector indicating wether a SNP have been used in the analysis (1's) or discarded (0's)}
  \item{S}{The joint emission probabilies for the current SNP and the privious SNP}
  \item{choose}{vector indivicating which privious SNP have been used to condition on}
  \item{mea}{The LD meassure for the <back> number of privious SNPs}
  \item{S1}{The emission probabilies for each single SNP}
  \item{hap}{The haplotypes probabilies for each privious SNPs}
  \item{maf}{the minor allele frequency for the SNPs}
	\item{path}{The famous viterbi path}
}
\references{ ~put references to the literature/web site here ~ }
\author{Anders Albrechtsen }
\note{ ~~further notes~~ 


}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
set.seed(3) 
n<-100 #number of individuals
snp=1000 #number of SNPs
freq<-runif(snp,min=0.01,max=0.99) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)
data<-matrix(rbinom(n*snp,2,freq),ncol=snp)+1
s<- sim_chr(snp,freq=freq, min=0.5, max=0.95, k=c(0.25,0.5,0.25), a=0.026, number_per_cm=5 )#simulate chromosomes for a sib pair
data<-rbind(t(s$geno)+1,data)
t<-relate(data,pair=c(1,2),pos=s$pos/1e6)#estimate tracts of relatedness
print(t)
plot(t)
points(s$pos/1e6,rep(0.4,s$snp),col=3-s$state,pch="|")#plot the true states
text(mean(max(s$pos)/2/1e6),0.45,"true States")
legend(150,0.9,paste("IBD=",2:0),col=1:3,lty=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
