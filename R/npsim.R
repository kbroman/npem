# npsim
#' Simulate data for the normal-Poisson mixture model
#'
#' Simulate data for a single plate using the normal-Poisson mixture model.
#'
#'
#' @param ests Parameter values as a vector for the form
#' (\eqn{\lambda_1}{lambda[1]},\ldots{},\eqn{\lambda_n}{lambda[n]}), a, b,
#' \eqn{\sigma}{sigma}
#' @param n Lengths of the \eqn{\lambda}{lambda} groups
#'
#' @return \item{k}{The simulated numbers of responding cells per well, useful
#' for comparison with those estimated from [npem.em()].}
#' \item{y}{The simulated transformed scintillation counts, to be used as input
#' for [npem.em()].}
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso [npem.em()]
#'
#' @references Broman et al. (1996) Estimation of antigen-responsive T cell
#' frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132
#'
#' @keywords datagen
#' @export
#' @importFrom stats rnorm rpois
#'
#' @examples
#'   ests <- c(0.5,1.5,2.5,5,15,5)
#'   n <- c(24,24,22)
#'   dat <- npsim(ests,n)
#'   out <- npem.em(dat$y,ests=ests,n=n)
#'   jitter <- runif(length(dat$k),-0.1,0.1)
#'   plot(dat$k + jitter, out$k, xlab="True no. cells",
#'        ylab="Estimated no. cells", lwd=2)
#'   plot(dat$y,out$k,type="n",xlab="Response",ylab="Estimated no. cells")
#'   for(i in unique(dat$k))
#'     points(dat$y[dat$k==i],out$k[dat$k==i],col=i+1,lwd=2)
#'
npsim<-
function(ests, n = c(24,24,24,22))
{
  if(length(n) != length(ests)-3) { # check number of groups
    if(length(n)==1) {
      n <- rep(n,length(ests)-3)
    }
    else {
      stop("Length of ests doesn't conform to length of n\n");
    }
  }
  n.groups <- length(n)

  ests[ests < 1e-10] <- 1e-10

  k <- vector("list",n.groups)
  for(i in 1:n.groups)
    k[[i]] <- rpois(n[i],ests[i])
  k <- unlist(k)
  y <- rnorm(sum(n),ests[n.groups+1]+ests[n.groups+2]*k,ests[n.groups+3])

  list(k = k, y = y)
}
