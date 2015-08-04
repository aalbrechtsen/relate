\name{plot.relate}
\alias{plot.relate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R generic plot function used for plotting relate.* functions }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.relate(x, col = 1:3, lwd = 2, ylab = "probability", xlab = "position (Mb)", chr = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{col}{ ~~Describe \code{col} here~~ }
  \item{lwd}{ ~~Describe \code{lwd} here~~ }
  \item{ylab}{ ~~Describe \code{ylab} here~~ }
  \item{xlab}{ ~~Describe \code{xlab} here~~ }
  \item{chr}{ ~~Describe \code{chr} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function (x, col = 1:3, lwd = 2, ylab = "probability", xlab = "position (Mb)", 
    chr = NULL, ...) 
{
    post <- t(x$post)
    pos <- x$position
    if (is.null(x$chr)) {
        plot(pos[1:(x$snp - 1)], post[2:(x$snp), 1], type = "l", 
            ylim = c(0, 1), col = col[1], lwd = lwd, ylab = ylab, 
            xlab = xlab, ...)
        lines(pos[1:(x$snp - 1)], post[2:(x$snp), 2], col = col[2], 
            lwd = lwd)
        lines(pos[1:(x$snp - 1)], post[2:(x$snp), 3], col = col[3], 
            lwd = lwd)
    }
    else {
        if (!is.null(chr)) {
            pos <- pos[x$chr == chr]
            post <- post[x$chr == chr, ]
            plot(pos[-length(pos)], post[-1, 1], type = "l", 
                ylim = c(0, 1), col = col[1], lwd = lwd, ylab = ylab, 
                xlab = xlab, ...)
            lines(pos[-length(pos)], post[-1, 2], col = col[2], 
                lwd = lwd)
            lines(pos[-length(pos)], post[-1, 3], col = col[3], 
                lwd = lwd)
        }
        else {
            C <- names(table(x$chr))
            m <- c(0, cumsum(tapply(pos, x$chr, max)))
            pos2 <- rep(NA, length(pos))
            for (tal in 1:length(C)) pos2[x$chr == C[tal]] <- pos[x$chr == 
                C[tal]] + m[tal]
            plot(pos2[1:(x$snp - 1)], post[2:(x$snp), 1], type = "l", 
                ylim = c(0, 1), col = col[1], lwd = lwd, ylab = ylab, 
                xlab = xlab, ...)
            lines(pos2[1:(x$snp - 1)], post[2:(x$snp), 2], col = col[2], 
                lwd = lwd)
            lines(pos2[1:(x$snp - 1)], post[2:(x$snp), 3], col = col[3], 
                lwd = lwd)
            abline(v = m)
            for (tal in 1:length(C)) text(m[tal] + diff(m)[tal]/2, 
                0.5, C[tal], col = "gray")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

