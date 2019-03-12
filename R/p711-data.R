#' Data for a limiting dilution assay
#'
#' This is data for a limiting dilution assay for a single patient, composed of
#' 6 pairs of plates at 6 different cell concentrations.  The data is from the
#' same subject (at the same time) as [p713.2()].
#'
#'
#' @format The data is a list with three components: \tabular{ll}{
#' `counts` \tab A list of length 12, each component of which is a vector
#' giving the square-root-transformed scintillation counts for a single plate.
#' \cr `cells` \tab A vector of length 12, giving the estimated number of
#' cells per well for each of the 12 plates. \cr `n` \tab Vector of length
#' 4, giving the number of wells per group, which is the same for each plate. }
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu} \cr
#' <https://github.com/kbroman/npem>
#' @seealso [p713.2()], [p711()], [npem.em()]
#' @source Michael Tigges, Chiron Biocine
#' @keywords datasets
"p711"
