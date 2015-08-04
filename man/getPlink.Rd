\name{getPlink}
\alias{getPlink}
\title{Read a PLINK binary data file as a list of chromosomes, genotypes and positions}
\description{
The package PLINK saves genome-wide association data in groups of three
files, with the extensions \code{.bed}, \code{.bim}, and
\code{.fam}. This function reads these files and return a list of chromosomes, positions and genotypes.
}
\usage{
getPlink(bed, bim, fam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{The name of the
    file containing the packed binary SNP genotype data}
  \item{bim}{The file containing the SNP descriptions}
  \item{fam}{The file containing subject (and, possibly, family)
    identifiers. This is basically a tab-delimited "pedfile"}
}
\details{
  If the \code{bed} argument does not contain a filename with the file
  extension \code{.bed}, then this extension is appended to the
  argument. The remaining two arguments are optional; their default
  values are obtained by replacing the \code{.bed} filename extension by
  \code{.bim} and \code{.fam} respectively. The code is very much inspired by the
	snpMatrix function called read.plink().
}
\value{
	A list that has the chromosomes, genotypes and positions. The SNP names are stripped.
}
\references{ ~put references to the literature/web site here ~ }
\author{Anders Albrechtsen \email{albrect@binf.ku.dk}}
\seealso{\code{\link{getPed}}
}
\examples{

  path<-paste(find.package("Relate"),"/data/500",sep="")
  plinkData<-getPlink(path)


}

\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
