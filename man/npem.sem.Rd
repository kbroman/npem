\name{npem.sem}
\alias{npem.sem}
\title{Obtain standard errors for the estimates from npem.em}
\description{
  Uses the SEM algorithm to obtain estimated standard errors for the
  MLEs obtained after fitting the normal-Poisson mixture model to data
  on a cell proliferation assay.
}
\usage{
npem.sem(y, npem.em.out, cells=10^6, start=npem.em.out$ests * 1.05 + 0.001,
         n=c(24, 24, 24, 22), n.plates=1, use.order.constraint=TRUE,
         all.se=TRUE, tol=1e-6, maxk=20, prnt=0, do.var=TRUE, maxit=1000)
}
\arguments{
  \item{y}{Vector of transformed scintillation counts, in
    lexicographical order (plate by plate and group by group within a
    plate.)}
  \item{npem.em.out}{Output from the function \code{\link{npem.em}}.}
  \item{cells}{Number of cells per well.  The \eqn{\lambda}{lambda}'s
    will be rescaled to give response per \eqn{10^6} cells.  This may be
    either a single number (if all wells have the same number of cells,
    or \eqn{10^6} if one wishes the \eqn{\lambda}{lambda}'s to not be
    rescaled), a value for each plate (vector of length \code{n.plates},
    or a value for each well (a vector of the same length as \code{y}).}
  \item{start}{Starting estimates, some small distance away from the
    MLE.  A vector of the form (\eqn{\lambda}{lambda}'s, (a, b,
    \eqn{\sigma}{sigma})'s).}
  \item{n}{Vector giving the number of wells within each group.  This
    may have length either n.groups (if all plates have the same
    number of wells per group) or n.groups*n.plates.}
  \item{n.plates}{The number of plates in the data.}
  \item{use.order.constraint}{If TRUE, force the constraint
    \eqn{\lambda_0 \le \lambda_i}{lambda[0] <= lambda[i]} for all
    \eqn{i \ge 1}{i >= 1}; otherwise, no constraints are applied.}
  \item{all.se}{If TRUE, do the full SEM algorithm; if FALSE, ignore the
    plate-specific parameters (a,b,\eqn{\sigma}{sigma})'s, to get
    estimated SEs for the \eqn{\lambda}{lambda}'s.}
  \item{tol}{Tolerance to determine when to stop the EM algorithm.}
  \item{maxk}{Maximum k value in sum calculating \eqn{E(k | y)}.}
  \item{prnt}{If 0, don't print anything; if 1, print out (at each step)
    the rate matrix and which elements have converged.}
  \item{do.var}{If TRUE, calculate the variance-covariance matrix and
    standard errors; if FALSE, only calculate the full-data information
    matrix and rate matrix.}
  \item{maxit}{Maximum number of iterations to perform.}
}
\details{
  Calculations are performed in a C routine.  It is important to first
  run \code{\link{npem.em}} with a very small value for \code{tol}, such
  as \eqn{10^{-13}}{10^-13}.
}
\value{
  \item{infor}{The full-data information matrix}
  \item{rates}{The rate matrix ("DM" in Meng and Rubin's notation).}
  \item{n.iter}{Number of iterations performed in calculating the rate
    matrix.}
  \item{var}{The estimated variance-covariance matrix.}
  \item{se}{The estimated standard errors.  (The square root of the
    diagnol of \code{var}.)}
}
\references{Broman et al. (1996) Estimation of antigen-responsive T cell
  frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132 \cr
  Dempster et al. (1977) Maximum likelihood estimation from
  incomplete data via the EM algorithm.  J Roy Statist Soc Ser B
  39:1-38 \cr
  Meng and Rubin (1991) Using EM to obtain asymptotic
  variance-covariance matrices: the SEM algorithm.  J Am Statist Asso
  86:899-909}
\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{https://github.com/kbroman/npem}}

\seealso{\code{\link{npem.em}}}

\examples{
  data(p713)
  start.pl3 <- npem.start(p713$counts[[3]],n=p713$n)
  out.pl3 <- npem.em(p713$counts[[3]],start.pl3,n=p713$n,tol=1e-13)
  out.sem.pl3 <- npem.sem(p713$counts[[3]],out.pl3,n=p713$n)
}

\keyword{models}
