\name{npsim}
\alias{npsim}
\title{Simulate data for the normal-Poisson mixture model}
\description{
  Simulate data for a single plate using the normal-Poisson mixture
  model.
}

\usage{
npsim(ests, n=c(24, 24, 24, 22))
}

\arguments{
  \item{ests}{Parameter values as a vector for the form
    (\eqn{\lambda_1}{lambda[1]},\ldots,\eqn{\lambda_n}{%
      lambda[n]}), a, b, \eqn{\sigma}{sigma}}
  \item{n}{Lengths of the \eqn{\lambda}{lambda} groups}
}
\value{
  \item{k}{The simulated numbers of responding cells per well, useful
    for comparison with those estimated from \code{\link{npem.em}}.}
  \item{y}{The simulated transformed scintillation counts, to be used as
  input for \code{\link{npem.em}}.}
}
\references{Broman et al. (1996) Estimation of antigen-responsive T cell
  frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132}
\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{https://github.com/kbroman/npem} }

\seealso{\code{\link{npem.em}}}

\examples{
  ests <- c(0.5,1.5,2.5,5,15,5)
  n <- c(24,24,22)
  dat <- npsim(ests,n)
  out <- npem.em(dat$y,ests=ests,n=n)
  jitter <- runif(length(dat$k),-0.1,0.1)
  plot(dat$k + jitter, out$k, xlab="True no. cells",
       ylab="Estimated no. cells", lwd=2)
  plot(dat$y,out$k,type="n",xlab="Response",ylab="Estimated no. cells")
  for(i in unique(dat$k))
    points(dat$y[dat$k==i],out$k[dat$k==i],col=i+1,lwd=2)
}

\keyword{datagen}
