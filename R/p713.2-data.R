#' Data for a limiting dilution assay
#'
#' This is data for a limiting dilution assay for a single patient, composed of
#' 5 pairs of plates at 5 different cell concentrations.  The data is from the
#' same subject (at the same time) as [p713()].
#'
#'
#' @format The data is a list with three components: \tabular{ll}{
#' `counts` \tab A list of length 10, each component of which is a vector
#' giving the square-root-transformed scintillation counts for a single plate.
#' \cr `cells` \tab A vector of length 10, giving the estimated number of
#' cells per well for each of the 10 plates. \cr `n` \tab Vector of length
#' 4, giving the number of wells per group, which is the same for each plate. }
#' @author Karl W Broman, \email{broman@@wisc.edu} \cr
#' <https://github.com/kbroman/npem>
#' @seealso [p713()], [p711()], [npem.em()]
#' @source Michael Tigges, Chiron Biocine
#' @keywords datasets
"p713.2"
