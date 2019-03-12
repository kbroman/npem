# npem.sem
#' Obtain standard errors for the estimates from npem.em
#'
#' Uses the SEM algorithm to obtain estimated standard errors for the MLEs
#' obtained after fitting the normal-Poisson mixture model to data on a cell
#' proliferation assay.
#'
#' Calculations are performed in a C routine.  It is important to first run
#' \code{\link{npem.em}} with a very small value for \code{tol}, such as
#' \eqn{10^{-13}}{10^-13}.
#'
#' @param y Vector of transformed scintillation counts, in lexicographical
#' order (plate by plate and group by group within a plate.)
#' @param npem.em.out Output from the function \code{\link{npem.em}}.
#' @param cells Number of cells per well.  The \eqn{\lambda}{lambda}'s will be
#' rescaled to give response per \eqn{10^6} cells.  This may be either a single
#' number (if all wells have the same number of cells, or \eqn{10^6} if one
#' wishes the \eqn{\lambda}{lambda}'s to not be rescaled), a value for each
#' plate (vector of length \code{n.plates}, or a value for each well (a vector
#' of the same length as \code{y}).
#' @param start Starting estimates, some small distance away from the MLE.  A
#' vector of the form (\eqn{\lambda}{lambda}'s, (a, b, \eqn{\sigma}{sigma})'s).
#' @param n Vector giving the number of wells within each group.  This may have
#' length either n.groups (if all plates have the same number of wells per
#' group) or n.groups*n.plates.
#' @param n.plates The number of plates in the data.
#' @param use.order.constraint If TRUE, force the constraint \eqn{\lambda_0 \le
#' \lambda_i}{lambda[0] <= lambda[i]} for all \eqn{i \ge 1}{i >= 1}; otherwise,
#' no constraints are applied.
#' @param all.se If TRUE, do the full SEM algorithm; if FALSE, ignore the
#' plate-specific parameters (a,b,\eqn{\sigma}{sigma})'s, to get estimated SEs
#' for the \eqn{\lambda}{lambda}'s.
#' @param tol Tolerance to determine when to stop the EM algorithm.
#' @param maxk Maximum k value in sum calculating \eqn{E(k | y)}.
#' @param prnt If 0, don't print anything; if 1, print out (at each step) the
#' rate matrix and which elements have converged.
#' @param do.var If TRUE, calculate the variance-covariance matrix and standard
#' errors; if FALSE, only calculate the full-data information matrix and rate
#' matrix.
#' @param maxit Maximum number of iterations to perform.
#'
#' @return \item{infor}{The full-data information matrix} \item{rates}{The rate
#' matrix ("DM" in Meng and Rubin's notation).} \item{n.iter}{Number of
#' iterations performed in calculating the rate matrix.} \item{var}{The
#' estimated variance-covariance matrix.} \item{se}{The estimated standard
#' errors.  (The square root of the diagnol of \code{var}.)}
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso \code{\link{npem.em}}
#'
#' @references Broman et al. (1996) Estimation of antigen-responsive T cell
#' frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132 \cr
#' Dempster et al. (1977) Maximum likelihood estimation from incomplete data
#' via the EM algorithm.  J Roy Statist Soc Ser B 39:1-38 \cr Meng and Rubin
#' (1991) Using EM to obtain asymptotic variance-covariance matrices: the SEM
#' algorithm.  J Am Statist Asso 86:899-909
#'
#' @keywords models
#' @export
#' @useDynLib npem, .registration=TRUE
#'
#' @examples
#'   data(p713)
#'   start.pl3 <- npem.start(p713$counts[[3]],n=p713$n)
#'   out.pl3 <- npem.em(p713$counts[[3]],start.pl3,n=p713$n,tol=1e-13)
#'   out.sem.pl3 <- npem.sem(p713$counts[[3]],out.pl3,n=p713$n)
#'
npem.sem<-
function(y, npem.em.out, cells=10^6, start=npem.em.out$ests * 1.05 + 0.001,
         n=c(24, 24, 24, 22), n.plates=1, use.order.constraint=TRUE,
         all.se=TRUE, tol=1e-6, maxk=20, prnt=0, do.var=TRUE, maxit=1000)
{

  n.groups <- length(start)-3*n.plates
  n.par <- n.groups + 3*n.plates
  if(all.se) n.se <- n.par
  else n.se <- n.groups

  # if plates are all the same, then n can be given as the
  # number of wells per group for a single plate, and must
  # then be expanded
  if(length(n) == n.groups) n <- rep(n,n.plates)

  if(length(y) != sum(n)) {
    cat("y is not the correct length\n")
    cat("length(y) =", length(y), "\tbut n gives", sum(n), "\n")
    return()
  }

  # if cells contains one cell density per plate, it is converted
  # to have one per well, where we assume y is given as plate 1,
  # plate 2, ...
  if(length(cells) != length(y))
    cells <- as.numeric(matrix(rep(cells,length(y)/length(cells)),
                 byrow=TRUE, nrow=length(y)/length(cells)))

#  if(!is.loaded("npem_sem")) {
#    lib.file <- file.path(paste("npem", .Platform$dynlib.ext, sep=""))
#    dyn.load(lib.file)
#    cat(paste(" -Loaded", lib.file), "\n")
#  }
#  if(!is.loaded(C.symbol("npem_sem"))) dyn.load("npem.so")

  z<-.C("npem_sem",
     as.double(y),
     as.integer(n.plates),
     as.double(cells/10^6),
     as.integer(n.groups),
     as.integer(n),
     as.integer(rep(0,length(n)+1)),
     as.integer(rep(0,n.plates)),
     as.integer(rep(0,n.groups)),
     as.integer(n.se),
     as.double(npem.em.out$e),
     as.double(start),
     as.double(npem.em.out$k),
     as.double(npem.em.out$ksq),
     infor=as.double(rep(0,(n.se)^2)),
     rates=as.double(rep(0,(n.se)^2)),
     var=as.double(rep(0,(n.se)^2)),
     se=as.double(rep(0,(n.se))),
     as.integer(as.numeric(do.var)),
     as.double(tol),
     as.integer(maxk),
     as.integer(prnt),
     as.integer(rep(0,10*(n.par)^2+10*n.par)),
     as.double(rep(0,20*(n.par)+20*(n.par)^2)),
     as.integer(as.numeric(use.order.constraint)),
     n.iter=as.integer(as.numeric(maxit)),
     PACKAGE="npem")

  z$infor<-matrix(z$infor,n.se,n.se)
  z$rates<-matrix(z$rates,n.se,n.se)
  x <- names(npem.em.out$ests)[1:n.se]
  if(do.var) {
    z$var <- matrix(z$var,n.se,n.se)
    dimnames(z$var) <- list(x,x)
    names(z$se) <- x
    return(list(infor=z$infor,rates=z$rates,n.iter=z$n.iter,var=z$var,se=z$se))
  }
  else return(list(infor=z$infor,rates=z$rates,n.iter=z$n.iter))
}
