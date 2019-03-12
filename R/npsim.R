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
