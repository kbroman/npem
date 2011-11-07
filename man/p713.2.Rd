\name{p713.2}
\alias{p713.2}

\title{Data for a limiting dilution assay}
\description{ 
  This is data for a limiting dilution assay for a single patient,
  composed of 5 pairs of plates at 5 different cell concentrations.
  The data is from the same subject (at the same time) as \code{\link{p713}}.
} 
\source{
  Michael Tigges, Chiron Biocine
}
\usage{data(p713.2)}
\format{
  The data is a list with three components:
  \tabular{ll}{
    \code{counts} \tab A list of length 10, each component of which is a
    vector giving the square-root-transformed scintillation counts for a
    single plate. \cr
    \code{cells} \tab A vector of length 10, giving the estimated number
    of cells per well for each of the 10 plates. \cr
    \code{n} \tab Vector of length 4, giving the number of wells per
    group, which is the same for each plate. }
}
\author{Karl W Broman,
\email{kbroman@biostat.wisc.edu} \cr
\url{http://biosun01.biostat.biostat.wisc.edu/~kbroman/software/npem.html} }

\seealso{\code{\link{p713}}, \code{\link{p711}}, \code{\link{npem.em}}}

\keyword{datasets}
