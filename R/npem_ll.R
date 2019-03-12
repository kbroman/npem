# npem.ll
#' Calculate the log likelihood for the normal-Poisson model
#'
#' Calculates the log likelihood at a single point in the parameter space, for
#' the normal-Poisson mixture model with data on a cell proliferation assay.
#'
#' Calculations are performed in a C routine.
#'
#' @param y Vector of transformed scintillation counts, in lexicographical
#' order (plate by plate and group by group within a plate.)
#' @param ests Value of the parameters at which to calculate the log
#' likelihood, as a vector of length n.groups + 3*n.plates, of the form
#' (\eqn{\lambda}{lambda}'s, (a, b, \eqn{\sigma}{sigma})'s), where
#' \eqn{\lambda}{lambda} is the average number of responding cells per
#' \eqn{10^6} cells for a group, and (a, b, \eqn{\sigma}{sigma}) are the
#' plate-specific parameters.
#' @param cells Number of cells per well.  The \eqn{\lambda}{lambda}'s will be
#' rescaled to give response per \eqn{10^6} cells.  This may be either a single
#' number (if all wells have the same number of cells, or \eqn{10^6} if one
#' wishes the \eqn{\lambda}{lambda}'s to not be rescaled), a value for each
#' plate (vector of length \code{n.plates}, or a value for each well (a vector
#' of the same length as \code{y}).
#' @param n Vector giving the number of wells within each group.  This may have
#' length either n.groups (if all plates have the same number of wells per
#' group) or n.groups*n.plates.
#' @param n.plates The number of plates in the data.
#' @param maxk Maximum k value in sum calculating \eqn{E(k | y)}.
#'
#' @return \item{loglik}{The log likelihood function calculated at the point
#' \code{ests} in the parameter space.}
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso \code{\link{npem.em}}
#'
#' @references Broman et al. (1996) Estimation of antigen-responsive T cell
#' frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132
#'
#' @keywords models
#' @export
#' @useDynLib npem, .registration=TRUE
#'
#' @examples
#'   data(p713)
#'   start.pl3 <- npem.start(p713$counts[[3]],n=p713$n)
#'   out.pl3 <- npem.em(p713$counts[[3]],start.pl3,n=p713$n)
#'   npem.ll(p713$counts[[3]],start.pl3,n=p713$n)
#'   npem.ll(p713$counts[[3]],out.pl3$ests,n=p713$n)
#'
npem.ll <-
function(y, ests, cells=10^6, n=c(24, 24, 24, 22), n.plates=1, maxk=30)
{

  n.groups <- length(ests) - 3*n.plates

  # if plates are all the same, then n can be given as the
  # number of wells per group for a single plate, and must
  # then be expanded
  if(length(n) == n.groups) n <- rep(n,n.plates)

  if(length(y) != sum(n)) {
    cat("y is not the correct length\n")
    cat("length(y) =", length(y), "\tbut n gives", sum(n), "\n")
    return()
  }

  len.group <- rep(0,n.groups)
  len.plate <- rep(0,n.plates)
  for(i in 1:length(len.group))
    len.group[i] <- sum(n[(i-1)+seq(1,length(n),by=n.groups)])
  for(i in 1:length(len.plate))
    len.plate[i] <- sum(n[(i-1)*n.groups+1:n.groups])

  # if cells contains one cell density per plate, it is converted
  # to have one per well, where we assume y is given as plate 1,
  # plate 2, ...
  if(length(cells) != length(y))
    cells <- as.numeric(matrix(rep(cells,length(y)/length(cells)),
                 byrow=TRUE, nrow=length(y)/length(cells)))

#  if(!is.loaded("npem_ll")) {
#    lib.file <- file.path(paste("npem_em", .Platform$dynlib.ext, sep=""))
#    dyn.load(lib.file)
#    cat(paste(" -Loaded", lib.file), "\n")
#  }
#  if(!is.loaded(C.symbol("npem_ll"))) dyn.load("npem_em.so")

  .C("npem_ll",
     as.double(y),
     as.integer(n.plates),
     as.double(cells/10^6),
     as.integer(n.groups),
     as.integer(n),
     as.integer(c(0,cumsum(n))),
     as.integer(len.plate),
     as.integer(len.group),
     as.double(ests),
     loglik=as.double(0),
     as.integer(maxk),
     PACKAGE="npem")$loglik

}
