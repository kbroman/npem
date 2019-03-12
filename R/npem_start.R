# npem.start
#' Obtain approximate starting values for npem.em
#'
#' Obtains crude starting values, needed for the EM algorithm.
#'
#' [I should describe the algorithm in more detail here.]
#'
#' @param y Vector of transformed scintillation counts, in lexicographical
#' order (plate by plate and group by group within a plate.)
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
#' @param n.groups The number of groups. (This is needed here but not
#' elsewhere, because usually I figure it out from `n.plates` and the
#' length of the argument `ests`.)
#' @param n.sd Number of SDs above the mean to use as a cutoff
#' @param cv Coefficient of variation (= SD/ave) used in randomizing the
#' starting point; use cv=0 to avoid randomization.
#'
#' @return \item{ests}{The parameter estimates to use as starting values for
#' the EM algorithm, as a vector of length n.groups + 3*n.plates, of the form
#' (\eqn{\lambda}{lambda}'s, (a, b, \eqn{\sigma}{sigma})'s), where
#' \eqn{\lambda}{lambda} is the average number of responding cells per
#' \eqn{10^6} cells for a group, and (a, b, \eqn{\sigma}{sigma}) are the
#' plate-specific parameters.}
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso [npem.em()]
#'
#' @references Broman et al. (1996) Estimation of antigen-responsive T cell
#' frequencies in PBMC from human subjects.  J Immunol Meth 198:119-132
#'
#' @keywords utilities
#' @export
#' @useDynLib npem, .registration=TRUE
#'
#' @examples
#'   data(p713)
#'   start.pl3 <- npem.start(p713$counts[[3]],n=p713$n)
#'   out.pl3 <- npem.em(p713$counts[[3]],start.pl3,n=p713$n,tol=1e-13)
#'
npem.start <-
function(y, cells=10^6, n=c(24, 24, 24, 22), n.plates=1, n.groups=4,
         n.sd = 2, cv=0.33)
{

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

  len.work <- n.plates*(2*n.groups+2)

#  if(!is.loaded("npem_start")) {
#    lib.file <- file.path(paste("npem", .Platform$dynlib.ext, sep=""))
#    dyn.load(lib.file)
#    cat(paste(" -Loaded", lib.file), "\n")
#  }
#  if(!is.loaded(C.symbol("npem_start"))) dyn.load("npem.so")

  z <- .C("npem_start",
     y = as.double(y),
     as.integer(n.plates),
     as.double(cells/10^6),
     as.integer(n.groups),
     as.integer(n),
     as.integer(rep(0,length(n)+1)),
     as.integer(rep(0,n.plates)),
     as.double(cv),
     as.double(n.sd),
     ests=as.double(rep(0,n.groups+3*n.plates)),
     as.double(rep(0,len.work)),
          PACKAGE="npem")

  z$ests
}
