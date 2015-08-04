\name{ld.snp3}
\alias{ld.snp3}

\title{ Function to calculate pairwise D', r-squared, but also haplotype
  frequencies}
\description{
  \code{ld.snp3} takes an genotype matrix and a suitable depth LOD and returns a list of
  r-square, lods and pba's.
}
\usage{
ld.snp3(snpdata, back = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpdata}{An genotype matrix with M samples of
    N snps}
  \item{back}{The depth or lag of pair-wise calculation. Should be
      between 1 and N-1; default 100.}
 }
\details{
	This function was used during development and is quite borrowed
  from the snpMatrix function called ld.snp().  The cubic equation and
  quadratic equation solver code is borrowed from GSL (GNU Scientific
  Library).
  
}
\value{
  \item{dprime}{D'}
  \item{rsq2}{$r^2$}
  \item{r}{signed $r$}
  \item{lod}{Log of Odd's}
  All the matrices are defined such that the ($n, m$)th entry is the
  pair-wise value between the ($n$)th snp and the $(n+m)$th snp. Hence
  the lower right triangles are always filled with zeros. (See example
  section for the actual layout)
  
  Invalid values are represented by an out-of-range value - currently we
  use -1 for D', $r^2$ (both of which are between 0 and 1), and -2 for
  $r$ (valid values are between -1 and +1). lod is set to zero in
  most of these invalid cases. (lod can be any value so it is not indicative).

}
\references{
  Anders albrechtsen, Thorfinn Sand Korneliussen
    Clayton, D.G. and Leung, Hin-Tak (2007) An R package for analysis of
  whole-genome association studies. 
  \emph{Human Heredity} \bold{64}:45-51.\cr
  GSL (GNU Scientific Library) \url{http://www.gnu.org/software/gsl/}
}
\author{Thorfinn Sand \email{thorfinn@binf.ku.dk}}
\keyword{dplot}
\keyword{htest}
\keyword{models}
