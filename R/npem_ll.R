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
