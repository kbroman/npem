
######################################################################
#
# npem
# version 0.50
# 22 August 2000
#
# Karl W Broman
#
# npem_em.R : npem.em, npem.ll, npsim
#
# July, 1995 (revised 10/11/95, 2/3/96, 8/21/2000, 6/22/2003)
#
######################################################################
######################################################################
#
#  npem.em
#
#  This function is the main S function to perform the EM algorithm
#  for the normal/Poisson model for a full NPEM.  It calls the C
#  function npem_em.  See the file "npem_em.h" for further details on
#  the model.
#
#       y = transformed scintillation counts, as a vector
#
#       ests = initial estimates (lambda's, (a, b, sigma)'s)
#
#       cells = cells/well, either of length "n.plates" or
#               of the same length as y, or of length 1
#
#       n = (n_i) = number of y's in the lambda groups
#          or = (n_{pi}) = number of y's in the lambda groups
#                          for all plates
#
#       n.plates = number of plates
#
#       use.order.constraint = T if you want lambda_1 <= lambda_i
#
#       maxit = maximum number of iterations
#
#       tol = tolerance value to determine when to stop the algorithm
#
#       maxk = maximum k value in sum in calculating E(k_ij|y_ij) in
#              the E-step
#
#       prnt = 0 : don't print anything out
#            = 1 : print the log likelihood at each iteration
#            = 2 : print the est's and the log lik. at each iteration
#
######################################################################

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

######################################################################
#
# npem.ll
#
# get log likelihood
#
######################################################################

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


######################################################################
#
#  npsim
#
#  This function simulates data for a single plate
#  using the Normal/Poisson model (described with the program "npem").
#
######################################################################

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


