# npem
#' Fit the normal-Poisson model to data on a cell proliferation assay
#'
#' Uses a version of the EM algorithm to fit the normal-Poisson mixture model
#' to data on a cell proliferation assay.
#'
#' Calculations are performed in a C routine.  \[I should describe the
#' normal-Poisson mixture model here.\]
#'
#' @param y Vector of transformed scintillation counts, in lexicographical
#' order (plate by plate and group by group within a plate.)
#' @param ests Initial parameter estimates, as a vector of length n.groups +
#' 3*n.plates, of the form (\eqn{\lambda}{lambda}'s, (a, b,
#' \eqn{\sigma}{sigma})'s), where \eqn{\lambda}{lambda} is the average number
#' of responding cells per \eqn{10^6} cells for a group, and (a, b,
#' \eqn{\sigma}{sigma}) are the plate-specific parameters.
#' @param cells Number of cells per well.  The \eqn{\lambda}{lambda}'s will be
#' rescaled to give response per \eqn{10^6} cells.  This may be either a single
#' number (if all wells have the same number of cells, or \eqn{10^6} if one
#' wishes the \eqn{\lambda}{lambda}'s to not be rescaled), a value for each
#' plate (vector of length `n.plates`, or a value for each well (a vector
#' of the same length as `y`).
#' @param n Vector giving the number of wells within each group.  This may have
#' length either n.groups (if all plates have the same number of wells per
#' group) or n.groups*n.plates.
#' @param n.plates The number of plates in the data.
#' @param use.order.constraint If TRUE, force the constraint \eqn{\lambda_0 \le
#' \lambda_i}{lambda[0] <= lambda[i]} for all \eqn{i \ge 1}{i >= 1}; otherwise,
#' no constraints are applied.
#' @param maxit Maximum number of EM iterations to perform.
#' @param tol Tolerance to determine when to stop the EM algorithm.
#' @param maxk Maximum k value in sum calculating E(k | y).
#' @param prnt If 0, don't print anything; if 1, print the log likelihood at
#' each iteration; and if 2, print the est's and the log lik. at each
#' iteration.
#'
#' @return \item{ests}{The estimated parameters in same form as the input
#' argument `ests`.} \item{k}{The estimated number of responding cells per
#' well, \eqn{E(k|y)}.} \item{ksq}{\eqn{E(k^2|y)}} \item{loglik}{The value of
#' the log likelihood at `ests`.} \item{n.iter}{Number of iterations
#' performed.}
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso [npem.sem()], [npem.start()],
#' [npsim()], [npem.ll()]
#'
#' @references Broman et al. (1996) Estimation of antigen-responsive T cell
#' frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132 \cr
#' Dempster et al. (1977) Maximum likelihood estimation from incomplete data
#' via the EM algorithm.  J Roy Statist Soc Ser B 39:1-38
#'
#' @keywords models
#' @export
#' @useDynLib npem, .registration=TRUE
#'
#' @examples
#'   # get access to an example data set
#'   data(p713)
#'
#'   # analysis of plate3
#'     # get starting values
#'   start.pl3 <- npem.start(p713$counts[[3]],n=p713$n)
#'     # get estimates
#'   out.pl3 <- npem.em(p713$counts[[3]],start.pl3,n=p713$n)
#'     # look at log likelihood at starting and ending points
#'   npem.ll(p713$counts[[3]],start.pl3,n=p713$n)
#'   npem.ll(p713$counts[[3]],out.pl3$ests,n=p713$n)
#'     # repeat with great precision, starting at previous endpoint
#'   out.pl3 <- npem.em(p713$counts[[3]],out.pl3$ests,
#'                      n=p713$n,tol=1e-13)
#'     # run SEM algorithm to get standard errors
#'   out.sem.pl3 <- npem.sem(p713$counts[[3]],out.pl3,n=p713$n)
#'   round(out.pl3$ests,3)
#'   round(out.sem.pl3$se,3)
#'
#'   # repeat the above for the pair, plates 3 and 4
#'     # get starting values
#'   start.pl34 <- npem.start(unlist(p713$counts[3:4]),n=p713$n,n.plates=2)
#'     # get estimates
#'   out.pl34 <- npem.em(unlist(p713$counts[3:4]),start.pl34,n=p713$n,n.plates=2)
#'     # look at log likelihood at starting and ending points
#'   npem.ll(unlist(p713$counts[3:4]),start.pl34,n=p713$n,n.plates=2)
#'   npem.ll(unlist(p713$counts[3:4]),out.pl34$ests,n=p713$n,n.plates=2)
#'     # repeat with great precision, starting at previous endpoint
#'   out.pl34 <- npem.em(unlist(p713$counts[3:4]),out.pl34$ests,
#'                      n=p713$n,tol=1e-13,n.plates=2)
#'     # run SEM algorithm to get standard errors
#'   out.sem.pl34 <- npem.sem(unlist(p713$counts[3:4]),out.pl34,n=p713$n,n.plates=2)
#'   round(out.pl34$ests,3)
#'   round(out.sem.pl34$se,3)
#'
npem.em <-
function(y, ests, cells=10^6, n=c(24, 24, 24, 22), n.plates=1,
         use.order.constraint=TRUE, maxit=2000, tol=1e-6, maxk=20, prnt=0)
{

  n.groups <- length(ests)-3*n.plates

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

#  if(!is.loaded("npem_em")) {
#    lib.file <- file.path(paste("npem", .Platform$dynlib.ext, sep=""))
#    dyn.load(lib.file)
#    cat(paste(" -Loaded", lib.file), "\n")
#  }

#  if(!is.loaded(C.symbol("npem_em"))) dyn.load("npem.so")

  z<-  .C("npem_em",
          as.double(y),
          as.integer(n.plates),
          as.double(cells/10^6),
          as.integer(n.groups),
          as.integer(n),
          as.integer(rep(0,length(n)+1)),
          as.integer(rep(0,n.plates)),
          as.integer(rep(0,n.groups)),
          ests=as.double(ests),
          k=as.double(rep(0,length(y))),
          ksq=as.double(rep(0,length(y))),
          loglik=as.double(0),
          n.iter=as.integer(maxit),
          as.double(tol),
          as.integer(maxk),
          as.integer(prnt),
          as.integer(rep(0,n.groups)),
          as.double(rep(0,3*n.groups+8*n.plates)),
          as.integer(as.numeric(use.order.constraint)),
          PACKAGE="npem")

  if(n.groups==4)
        names(z$ests) <- c(paste("lam",c("ca","d","b","t"),sep="."),
                paste(c("a","b","sigma"),sort(rep(1:n.plates,3)),
                sep=""))
  else
        names(z$ests) <- c(paste("lam",1:n.groups,sep="."),
                paste(c("a","b","sigma"),sort(rep(1:n.plates,3)),
                sep=""))

  list(ests=z$ests,k=z$k,ksq=z$ksq,loglik=z$loglik,n.iter=z$n.iter)

}
