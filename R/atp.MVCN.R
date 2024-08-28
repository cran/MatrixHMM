#' Atypical Detection Points Using Matrix-Variate Contaminated Normal Hidden Markov Models
#'
#' Detects atypical matrices via matrix-variate contaminated normal Hidden Markov Models.
#'
#' @param Y An array with dimensions \code{p} x \code{r} x \code{num} x \code{t}, where \code{p} is the number of
#'     variables in the rows of each data matrix, \code{r} is the number of variables in the columns of each
#'     data matrix, \code{num} is the number of data observations, and \code{t} is the number of time points.
#' @param pgood An array with dimensions \code{num} x \code{t} x \code{k} containing the estimated probability of being typical for each point, given the time and state.
#' @param class An \code{num} x \code{t} matrix containing the state memberships.
#'
#' @return An \code{num} x \code{t} matrix containing, for each observation and time, a 0 if it that matrix is typical and 1 otherwise.
#' @export
#' @examples
#' data("simData2")
#' Y <- simData2$Y

#' init <- Eigen.HMM_init(Y = Y, k = 2, density = "MVCN", mod.row = "EEE", mod.col = "EE", nstartR = 1)
#' fit <- Eigen.HMM_fit(Y = Y, init.par = init, nThreads = 1)
#' atp <- atp.MVCN(Y = Y,
#'                pgood = fit[["results"]][[1]][[1]][[1]][["pgood"]],
#'                class = fit[["results"]][[1]][[1]][[1]][["class"]])
#' which(atp==1)
#' which(simData2[["atp.tr"]]==1)
atp.MVCN <- function(Y, pgood, class) {
  num <- dim(Y)[3]
  t <- dim(Y)[4]

  w.obs <- matrix(0, nrow = num, ncol = t)

  for (j in 1:num) {
    for (l in 1:t) {
      w.obs[j, l] <- pgood[j, l, class[j, l]]
    }
  }

  innerCN <- matrix(0, nrow = num, ncol = t)

  for (j in 1:num) {
    for (l in 1:t) {
      innerCN[j, l] <- ifelse(w.obs[j, l] < 0.5, 1, 0) # 1 is atypical, 0 is typical
    }
  }

  return(atp = innerCN)
}
