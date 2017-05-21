######################################################################
#
# npem
# version 0.50
# 22 August 2000
#
# Karl W Broman
#
# npem_sem.R : npem.sem
#
# July, 1995 (revised 10/11/95, 2/3/96, 3/12/96, 4/30/96, 8/21/2000,
#                     6/22/2003)
#
######################################################################
######################################################################
#
#  npem.sem
#
#  This is the main S function which applies the SEM algorithm in
#  order to calculate an estimated variance-covariance matrix for the
#  normal/Poisson model with a full NPEM.  It calls the function
#  "npem_sem"  See the file "npem.h" for further details on the
#  model.
#
#  y, cells, n, n.plates use.order.constraint, tol, maxk
#       (similar to those for the program "npem.em")
#
#  npem.em.out = output from the program "npem.em"
#
#  start = starting parameter values for the SEM algorithm
#
#  all.se = T  -> get all SEs
#           F  -> ignore the (a,b,sigma)
#
#  do.var = T  -> calculate var-cov matrix
#           F  -> just calculate full-data information matrix and
#                 the rate matrix
#
#  prnt   = 0  -> don't print anything out
#         = 1  -> print out the rate matrix and which elements have
#                 converged at each step
#
######################################################################

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
