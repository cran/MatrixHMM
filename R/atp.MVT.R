#' Atypical Detection Points Using Matrix-Variate t Hidden Markov Models
#'
#' Detects atypical matrices via matrix-variate t Hidden Markov Models given a specified value of \code{epsilon}.
#'
#' @param Y An array with dimensions \code{p} x \code{r} x \code{num} x \code{t}, where \code{p} is the number of
#'     variables in the rows of each data matrix, \code{r} is the number of variables in the columns of each
#'     data matrix, \code{num} is the number of data observations, and \code{t} is the number of time points.
#' @param M An array with dimensions \code{p} x \code{p} x \code{k}, where \code{k} is the number of states, containing the mean matrices.
#' @param U An array with dimensions \code{p} x \code{p} x \code{k}, where \code{k} is the number of states, containing the row covariance (scale) matrices.
#' @param V An array with dimensions \code{r} x \code{r} x \code{k}, where \code{k} is the number of states, containing the column covariance (scale) matrices.
#' @param class An \code{num} x \code{t} matrix containing the state memberships.
#' @param epsilon A numeric value specifying the selected percentile of the chi-squared distribution with \code{pr} degrees of freedom.
#'
#' @return An \code{num} x \code{t} matrix containing, for each observation and time, a 0 if it that matrix is typical and 1 otherwise.
#' @export
#' @examples
#'data("simData2")
#'Y <- simData2$Y

#'init <- Eigen.HMM_init(Y = Y, k = 2, density = "MVT", mod.row = "EEE", mod.col = "EE", nstartR = 1)
#'fit <- Eigen.HMM_fit(Y = Y, init.par = init, nThreads = 1)
#'atp <- atp.MVT(Y = Y, M = fit[["results"]][[1]][[1]][[1]][["M"]],
#'               U = fit[["results"]][[1]][[1]][[1]][["U"]],
#'               V = fit[["results"]][[1]][[1]][[1]][["V"]],
#'               class = fit[["results"]][[1]][[1]][[1]][["class"]],
#'               epsilon = 0.99)
#'which(atp==1)
#'which(simData2[["atp.tr"]]==1)
atp.MVT <- function(Y, M, U, V, class, epsilon) {
  tr <- function(x) {
    return(sum(diag(x)))
  }

  p <- dim(Y)[1]
  r <- dim(Y)[2]
  num <- dim(Y)[3]
  t <- dim(Y)[4]
  k <- dim(M)[3]

  Yresh <- array(Y, dim = c(p, r, num * t))
  deltas <- matrix(NA, nrow = num * t, ncol = k)

  for (j in 1:k) {
    deltas[, j] <- sapply(1:(num * t), function(l) tr(solve(U[, , j]) %*% (Yresh[, , l] - M[, , j]) %*% solve(V[, , j]) %*% t(Yresh[, , l] - M[, , j])))
  }

  innerT <- numeric(num * t)

  for (j in 1:(num * t)) {
    innerT[j] <- ifelse((1 - stats::pchisq(deltas[j, class[j]], df = p * r)) < (1 - epsilon), 1, 0) # 1 is atypical, 0 is typical
  }

  innerT2 <- matrix(innerT, num, t)

  return(atp = innerT2)
}
