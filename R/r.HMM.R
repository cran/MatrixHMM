#' Random Number Generation for Matrix-Variate Hidden Markov Models
#'
#' Generates random numbers for matrix-variate Hidden Markov Models (HMMs) based on matrix-variate normal, t, and contaminated normal distributions.
#'
#' @param density A character string specifying the distribution to use for the HMM.
#'                Possible values are: "MVN" for the matrix-variate normal distribution, "MVT" for the matrix-variate t-distribution,
#'                and "MVCN" for the matrix-variate contaminated normal distribution.
#' @param num An integer specifying the number of random matrices to generate.
#' @param t An integer specifying the number of time points.
#' @param PI A matrix representing the transition probability matrix.
#' @param M An array with dimensions \code{p} x \code{r} x \code{k}, where \code{k} is the number of states, containing the mean matrices.
#' @param U An array with dimensions \code{p} x \code{p} x \code{k}, where \code{k} is the number of states, containing the row covariance (scale) matrices.
#' @param V An array with dimensions \code{r} x \code{r} x \code{k}, where \code{k} is the number of states, containing the column covariance (scale) matrices.
#' @param IP A numeric vector of length \code{k} containing the initial probability weights.
#' @param nu A numeric vector of length \code{k} containing the degrees of freedom for each state in the MVT distribution.
#' @param alpha A numeric vector of length \code{k} containing the proportion of typical points in each state for the MVCN distribution.
#' @param eta A numeric vector of length \code{k} containing the inflation parameters for each state in the MVCN distribution.
#'
#' @return A list containing the following elements:
#' \item{Y}{An array with dimensions \code{p} x \code{r} x \code{num} x \code{t} containing the generated data.}
#' \item{obs.states}{An \code{num} x \code{t} matrix containing the state memberships.}
#' @export
#' @examples
#' p <- 2
#' r <- 3
#' num <- 50
#' t <- 3
#' k <- 2
#' IP <- c(0.5, 0.5)
#' PI <- matrix(c(0.9, 0.1, 0.3, 0.7), nrow = k, ncol = k, byrow = TRUE)
#' M <- array(NA, dim = c(p, r, k))
#' M[,,1]<- matrix(c(0,1,1,
#'                  -1,-1.5,-1),nrow = p, ncol = r, byrow = TRUE)
#' M[,,2]<- M[,,1]+3
#' U <- array(NA, dim = c(p, p, k))
#' V <- array(NA, dim = c(r, r, k))
#' U[, , 1] <- U[, , 2] <- matrix(c(1.73, -0.59, -0.59, 2.52), nrow = p, ncol = p, byrow = TRUE)
#' V[, , 1] <- V[, , 2] <- matrix(c(0.69, 0.23, -0.03,
#'                                  0.23, 0.48,  0.16,
#'                                 -0.03, 0.16,  0.88), nrow = r, ncol = r, byrow = TRUE)
#' nu <- c(4.5, 6.5)
#' simData <- r.HMM(density = "MVT", num = num, t = t, PI = PI,
#'                 M = M, U = U, V = V, IP = IP, nu = nu)
r.HMM <- function(density, num, t, PI, M, U, V, IP, nu, alpha, eta) {
  run.mc.sim <- function(P, times, IP) {
    # number of possible states
    num.states <- nrow(P)

    # stores the states X_t through time
    states <- numeric(times)

    # initialize variable for first state

    k <- length(IP)

    states[1] <- sample(1:k, size = 1, prob = IP)

    for (t in 2:times) {
      # probability vector to simulate next state X_{t+1}
      p <- P[states[t - 1], ]

      ## draw from multinomial and determine state
      states[t] <- which(stats::rmultinom(1, 1, p) == 1)
    }
    return(states)
  }
  rTMV <- function(n, M, U, V, nu) {
    w <- stats::rgamma(n = n, nu / 2, nu / 2)
    X <- array(0, dim = c(nrow(M), ncol(M), n))

    for (i in 1:n) {
      X[, , i] <- LaplacesDemon::rmatrixnorm(M = M, U = U / w[i], V = V)
    }

    return(X)
  }

  row.Y <- dim(U)[1]
  col.Y <- col.X <- dim(V)[2]
  k <- nrow(PI)
  good <- matrix(NA, num, ncol = t)

  Y <- array(NA, dim = c(row.Y, col.Y, num, t))

  # each column stores the sequence of states for a single chains
  chain.states <- matrix(NA, ncol = num, nrow = t)

  # simulate chains
  for (c in seq_len(num)) {
    chain.states[, c] <- run.mc.sim(PI, t, IP = IP)
  }

  if (density == "MVN") {
    for (j in seq_len(t)) {
      for (i in 1:num) {
        Y[, , i, j] <- LaplacesDemon::rmatrixnorm(M = M[, , chain.states[j, i]], U = U[, , chain.states[j, i]], V = V[, , chain.states[j, i]])
      }
    }
  } else if (density == "MVT") {
    for (j in seq_len(t)) {
      for (i in 1:num) {
        Y[, , i, j] <- rTMV(n = 1, M = M[, , chain.states[j, i]], U = U[, , chain.states[j, i]], V = V[, , chain.states[j, i]], nu = nu[chain.states[j, i]])
      }
    }
  } else if (density == "MVCN") {
    for (j in seq_len(t)) {
      for (i in 1:num) {
        good[i, j] <- stats::rbinom(n = 1, size = 1, prob = alpha[chain.states[j, i]])
      }
    }

    for (j in seq_len(t)) {
      for (i in 1:num) {
        if (good[i, j] == 1) {
          Y[, , i, j] <- LaplacesDemon::rmatrixnorm(M = M[, , chain.states[j, i]], U = U[, , chain.states[j, i]], V = V[, , chain.states[j, i]])
        } else {
          Y[, , i, j] <- LaplacesDemon::rmatrixnorm(M = M[, , chain.states[j, i]], U = eta[chain.states[j, i]] * U[, , chain.states[j, i]], V = V[, , chain.states[j, i]])
        }
      }
    }
  }

  chain.states <- t(chain.states)
  colnames(chain.states) <- paste("time", 1:t, sep = " ")

  if (density == "MVCN") {
    return(list(Y = Y, obs.states = chain.states, pgood = good))
  } else {
    return(list(Y = Y, obs.states = chain.states))
  }
}
