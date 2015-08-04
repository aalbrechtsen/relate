\name{sim_chr}
\alias{sim_chr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Program used to generate data }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_chr(snp, freq = 0.2, min = 0.5, max = 0.95, k = k, a = 3.854, number_per_cm = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ ~~Describe \code{snp} here~~ }
  \item{freq}{ ~~Describe \code{freq} here~~ }
  \item{min}{ ~~Describe \code{min} here~~ }
  \item{max}{ ~~Describe \code{max} here~~ }
  \item{k}{ ~~Describe \code{k} here~~ }
  \item{a}{ ~~Describe \code{a} here~~ }
  \item{number_per_cm}{ ~~Describe \code{number_per_cm} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 


}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (snp, freq = 0.2, min = 0.5, max = 0.95, k = k, a = 3.854, 
    number_per_cm = 1) 
{
    hap11 <- unrelated(snp, freq, min, max)
    hap12 <- unrelated(snp, freq, min, max)
    hap21 <- unrelated(snp, freq, min, max)
    hap22 <- unrelated(snp, freq, min, max)
    pos <- as.integer(rnorm(snp, mean = CM, sd = CM/10))
    pos <- round(pos/number_per_cm, 0)
    pos <- cumsum(pos)
    state = recom(pos/CM, k = k, a = a)
    hap11[state == 2] <- hap21[state == 2]
    hap12[state == 2] <- hap22[state == 2]
    hap11[state == 1] <- hap21[state == 1]
    geno <- cbind(hap11 + hap12, hap21 + hap22)
    obs <- list(geno = geno, state = state, k = k, a = a, freq = freq, 
        number_per_cm = number_per_cm, min = min, max = max, 
        snp = snp, pos = pos)
    class(obs) = "sim_chr"
    return(obs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

