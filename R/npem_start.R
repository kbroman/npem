
######################################################################
#
# npem
# version 0.50
# 22 August 2000
#
# Karl W Broman
#
# npem_start.R: npem.start
#
# 3/13/96 (revised 4/25/96, 8/22/2000, 6/22/2003)
#
######################################################################
######################################################################
#
# npem.start
#
# This function calls a C function, npem_start, in order to get a
# crude starting point for the EM algorithm.
#
# Input:
#
#   y = transformed scintillation counts
#
#   cells = number of cells per well (generally is given as 10^6)
#
#   n = (n_i) = number of y's in the lambda groups
#
#   n.plates = number of plates
#
#   n.sd = number of SDs above the mean to use as cutoff
#
#   cv = common coefficient of variation (= SD/ave) used in randomizing
#        starting points; use cv=0 to not randomize
#
######################################################################

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
