\name{p711}
\alias{p711}

\title{Data for a limiting dilution assay}
\description{
  This is data for a limiting dilution assay for a single patient,
  composed of 6 pairs of plates at 6 different cell concentrations.
}
\source{
  Michael Tigges, Chiron Biocine
}
\usage{data(p711)}
\format{
  The data is a list with three components:
  \tabular{ll}{
    \code{counts} \tab A list of length 12, each component of which is a
    vector giving the square-root-transformed scintillation counts for a
    single plate. \cr
    \code{cells} \tab A vector of length 12, giving the estimated number
    of cells per well for each of the 12 plates. \cr
    \code{n} \tab Vector of length 4, giving the number of wells per
    group, which is the same for each plate. }
}
\author{Karl W Broman,
\email{kbroman@biostat.wisc.edu} \cr
\url{https://github.com/kbroman/npem} }

\seealso{\code{\link{p713}}, \code{\link{p713.2}}, \code{\link{npem.em}}}

\keyword{datasets}
