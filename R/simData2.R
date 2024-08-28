#' A Simulated Dataset with Atypical Matrices
#'
#' A simulated dataset containing atypical matrices.
#' The data are initially generated from a matrix-variate normal Hidden Markov Model with 2 states and an EE - EE covariance structure.
#' Atypical matrices are then introduced by randomly replacing some of the original matrices with values from a uniform distribution.
#'
#' @usage data(simData2)
#' @format A list containing three elements:
#' \describe{
#'   \item{1)}{An array with \code{p = 2} variables in the rows, \code{r = 3} variables in the columns, \code{num = 50} matrices, and \code{t = 3} time points.}
#'   \item{2)}{An \code{num} x \code{t} matrix containing the state memberships.}
#'   \item{3)}{An \code{num} x \code{t} matrix identifying the atypical matrices, where atypical matrices are coded with a 1.}
#' }
"simData2"
