\name{getPed}
\alias{getPed}

\title{A function to extract values from a snpMatrix object,snp.support
  and subject.support }
\description{
  This functions main purpose is only for debugging purposes.
  It solely returns the genotypes,positions and chromosomes of a
  snpMatrix pedfile object. Which in turn consist of a list with 3 elements.
  }
  \usage{
#library(snpMatrix)
#ped <- read.snps.pedfile("data.ped")
vars <- getPed(pedfile)
}
\arguments{
  \item{pedfile}{Return value of read.snps.pedfile of snpMatrix }
}

\value{
  \item{chromosome}{This is an Integer containing the 'chromosome'}
  \item{position}{Vector of real values, containing the relative
    position to the chromosome}
  \item{data}{snp data genotypes}
}
\references{http://staff.pubhealth.ku.dk/~ande/web/software/relate.html }
\author{albrecht@binf.ku.dk, thorfinn@binf.ku.dk}

\seealso{ \code{\link{getPlink}}.}
\examples{
#library("snpMatrix")
##
##  SNP MATRIX MUST BE INSTALLED
#  path<-paste(.find.package("Relate"),"/data/500.ped",sep="")
#   ped<-read.snps.pedfile(path)
#  datas<-getPed(ped)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

