#' Initialization for ECM Algorithms in Matrix-Variate Hidden Markov Models
#'
#' Initializes the ECM algorithms used for fitting parsimonious matrix-variate Hidden Markov Models (HMMs).
#' Parallel computing is implemented and highly recommended for faster computations.
#'
#' @param Y An array with dimensions \code{p} x \code{r} x \code{num} x \code{t}, where \code{p} is the number of
#'     variables in the rows of each data matrix, \code{r} is the number of variables in the columns of each
#'     data matrix, \code{num} is the number of data observations, and \code{t} is the number of time points.
#' @param k An integer or vector indicating the number of states in the model(s).
#' @param density A character string specifying the distribution to use in the HMM.
#'                Possible values are: "MVN" for the matrix-variate normal distribution,
#'                "MVT" for the matrix-variate t-distribution, and "MVCN" for the matrix-variate contaminated normal distribution.
#' @param mod.row A character string indicating the parsimonious structure of the row covariance (or scale) matrices.
#'                Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV", or "all".
#'                When "all" is specified, all 14 parsimonious structures are considered.
#' @param mod.col A character string indicating the parsimonious structure of the column covariance (or scale) matrices.
#'                Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV", or "all".
#'                When "all" is specified, all 7 parsimonious structures are considered.
#' @param nstartR An integer specifying the number of random starts to consider.
#' @param nThreads A positive integer indicating the number of cores to use for parallel processing.
#' @param verbose A logical value indicating whether to display the running output.
#' @param seed A positive integer specifying the seed for random generation.
#'
#' @return A list containing the following elements:
#' \item{results}{A list of the results from the initialization.}
#' \item{k}{The number of states fitted in each model.}
#' \item{req.model}{A data frame listing the models that were initialized.}
#' \item{init.used}{A data frame listing the initializations used for the required models.}
#' \item{index}{A numeric vector to be used by the \code{Eigen.HMM_fit()} function.}
#' \item{dens}{The density used for the HMMs.}
#' @export
#' @importFrom foreach %dopar%
#' @examples
#' data(simData)
#' Y <- simData$Y
#' init <- Eigen.HMM_init(Y = Y, k = 2, density = "MVT", mod.row = "EEE", mod.col = "EE", nstartR = 10)
Eigen.HMM_init <- function(Y, k, density, mod.row = "all", mod.col = "all", nstartR = 50, nThreads = 1, verbose = FALSE, seed = 3) {
  r_Pars_init <- function(Y, k, density, mod.row, mod.col, nstartR = 50, seed = NULL) {
    dMVnorm <- function(X, M, U, V) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

      return(pdf)
    }
    dMVT <- function(X, M, U, V, nu) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      pdfvar <- (1 + delta / nu)^(-0.5 * ((p * r) + nu))
      pdfconst <- (det(U)^(-r / 2) * det(V)^(-p / 2) * gamma(0.5 * ((p * r) + nu))) / ((pi * nu)^(0.5 * (p * r)) * gamma(nu * 0.5))

      PDF <- pdfconst * pdfvar

      return(PDF)
    }
    dMVcont <- function(X, M, U, V, alpha, eta) {
      dMVnorm <- function(X, M, U, V) {
        num <- dim(X)[3] # sample size
        p <- nrow(X) # rows of X
        r <- ncol(X) # columns of X

        if (is.na(num)) {
          X <- as.matrix(X)
          delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
        } else {
          delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
        }

        pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

        return(pdf)
      }

      pdf <- alpha * dMVnorm(X = X, M = M, U = U, V = V) + (1 - alpha) * dMVnorm(X = X, M = M, U = eta * U, V = V)

      return(pdf)
    }

    tr <- function(x) {
      return(sum(diag(x)))
    }
    transit <- function(Y) {
      # Y is a matrix N \times T

      N <- dim(Y)[1]
      T <- dim(Y)[2]
      K <- max(Y)

      if (K == 1) {
        PI <- matrix(1, 1, 1)
      }

      if (K > 1) {
        PI <- matrix(0, nrow = K, ncol = K)

        for (i in 1:N) {
          for (t in 2:T) {
            PI[Y[i, t - 1], Y[i, t]] <- PI[Y[i, t - 1], Y[i, t]] + 1
          }
        }
        PI <- diag(1 / rowSums(PI)) %*% PI
      }

      return(PI)
    }

    # Dimensions

    num0 <- dim(Y)[3]
    p <- nrow(Y)
    r <- ncol(Y)
    t <- dim(Y)[4]

    num <- num0 * t
    Yresh <- array(Y, dim = c(p, r, num))

    # Create some objects

    prior <- matrix(NA, nstartR, k)
    llk <- rep(NA, nstartR)

    M <- array(0, dim = c(p, r, k, nstartR))

    WRY <- array(0, dim = c(p, p, k, nstartR))
    WCY <- array(0, dim = c(r, r, k, nstartR))

    sigmaUY <- array(0, dim = c(p, p, k, nstartR))
    sigmaVY <- array(0, dim = c(r, r, k, nstartR))

    sel <- array(0, dim = c(p, r, k))
    dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
    post2 <- array(NA, c(k, k, num0, t - 1)) # transition probabilities
    f <- array(NA, c(num0, t, k, nstartR))
    A <- B <- array(NA, c(num0, t, k, nstartR))
    PI <- array(NA, dim = c(k, k, nstartR))

    nu <- alpha <- eta <- rep(0, k)

    if (density == "MVT") {
      nu <- rep(30, k)
    } else if (density == "MVCN") {
      alpha <- rep(0.99, k)
      eta <- rep(1.01, k)
    }

    ## Random initialization ##

    eu <- matrix(0, nrow = num, ncol = k)
    classy <- numeric(num)
    rand.start <- matrix(0, nstartR, k)

    withr::with_seed(seed, for (i in 1:nstartR) {
        rand.start[i, ] <- sample(c(1:num), k)
     })

    evs <- numeric(nstartR)

    for (h in 1:nstartR) {
      skip_to_next <- FALSE

      ### part 0 ###

      tryCatch(
        {
          sec <- rand.start[h, ]

          for (j in 1:k) {
            sel[, , j] <- Yresh[, , sec[j]]
          }

          for (j in 1:k) {
            for (i in 1:(num)) {
              eu[i, j] <- norm((Yresh[, , i] - sel[, , j]), type = "F")
            }
          }

          for (i in 1:(num)) {
            classy[i] <- which.min(eu[i, ])
          }

          z <- mclust::unmap(classy)

          ### part 1 ###

          # M + Partial U & V #

          for (j in 1:k) {
            M[, , j, h] <- rowSums(Yresh * z[, j][slice.index(Yresh, 3)], dims = 2) / sum(z[, j])
            MX <- array(M[, , j, h], dim = c(p, r, num))

            WRY[, , j, h] <- tensor::tensor(aperm(tensor::tensor((Yresh - MX) * z[, j][slice.index(Yresh - MX, 3)], diag(r), 2, 1), c(1, 3, 2)), aperm((Yresh - MX), c(2, 1, 3)), c(2, 3), c(1, 3))
            WCY[, , j, h] <- tensor::tensor(aperm(tensor::tensor(aperm(Yresh - MX, c(2, 1, 3)) * z[, j][slice.index(aperm(Yresh - MX, c(2, 1, 3)), 3)], solve(WRY[, , j, h]), 2, 1), c(1, 3, 2)), (Yresh - MX), c(2, 3), c(1, 3))
          }

          # Row covariance matrix Y #

          if (mod.row == "EII") {
            if (k == 1) {
              phiY <- tr(WRY[, , , h]) / ((num) * p * r)
            } else {
              phiY <- tr(rowSums(WRY[, , , h], dims = 2)) / ((num) * p * r)
            }

            for (j in 1:k) {
              sigmaUY[, , j, h] <- phiY * diag(1, p, p)
            }
          }

          if (mod.row == "EEI") {
            if (k == 1) {
              deltaU <- diag(diag(WRY[, , , h]), p, p) / (det(diag(diag(WRY[, , , h]), p, p)))^(1 / p)

              phiY <- (det(diag(diag(WRY[, , , h]), p, p)))^(1 / p) / ((num) * r)
            } else {
              deltaU <- diag(diag(rowSums(WRY[, , , h], dims = 2)), p, p) / (det(diag(diag(rowSums(WRY[, , , h], dims = 2)), p, p)))^(1 / p)

              phiY <- (det(diag(diag(rowSums(WRY[, , , h], dims = 2)), p, p)))^(1 / p) / ((num) * r)
            }

            for (j in 1:k) {
              sigmaUY[, , j, h] <- phiY * deltaU
            }
          }

          if (mod.row == "EEE") {
            if (k == 1) {
              sigmaUY[, , 1, h] <- WRY[, , , h] / (num * r)
            } else {
              for (j in 1:k) {
                sigmaUY[, , j, h] <- rowSums(WRY[, , , h], dims = 2) / (num * r)
              }
            }
          }

          # Column covariance matrix Y #

          if (mod.col == "II") {
            for (j in 1:k) {
              sigmaVY[, , j, h] <- diag(1, r, r)
            }
          }

          if (mod.col == "EI") {
            if (k == 1) {
              deltaVY <- diag(diag(WCY[, , , h]), r, r) / (det(diag(diag(WCY[, , , h]), r, r)))^(1 / r)
            } else {
              deltaVY <- diag(diag(rowSums(WCY[, , , h], dims = 2)), r, r) / (det(diag(diag(rowSums(WCY[, , , h], dims = 2)), r, r)))^(1 / r)
            }

            for (j in 1:k) {
              sigmaVY[, , j, h] <- deltaVY
            }
          }

          if (mod.col == "EE") {
            if (k == 1) {
              sigmaVY[, , 1, h] <- WCY[, , , h] / (det(WCY[, , , h]))^(1 / r)
            } else {
              for (j in 1:k) {
                sigmaVY[, , j, h] <- rowSums(WCY[, , , h], dims = 2) / ((det(rowSums(WCY[, , , h], dims = 2)))^(1 / r))
              }
            }
          }

          # Density computation #

          if (density == "MVN") {
            for (j in 1:k) {
              dens[, j] <- dMVnorm(X = Yresh, M = M[, , j, h], U = sigmaUY[, , j, h], V = sigmaVY[, , j, h])
            }
          } else if (density == "MVT") {
            for (j in 1:k) {
              dens[, j] <- dMVT(X = Yresh, M = M[, , j, h], U = sigmaUY[, , j, h], V = sigmaVY[, , j, h], nu = nu[j])
            }
          } else if (density == "MVCN") {
            for (j in 1:k) {
              dens[, j] <- dMVcont(X = Yresh, M = M[, , j, h], U = sigmaUY[, , j, h], V = sigmaVY[, , j, h], alpha = alpha[j], eta = eta[j])
            }
          }

          if (k == 1) {
            prior[h, ] <- 1
          } else {
            prior[h, ] <- colMeans(z)
          }

          ### part 2 - in log appraoch to prevent underflow ###

          f[, , , h] <- array(dens, c(num0, t, k))

          PI[, , h] <- transit(matrix(classy, num0, t))

          foo <- matrix(rep(prior[h, ], each = num0), ncol = k) * f[, 1, , h]
          sumfoo <- rowSums(foo)
          lscale <- log(sumfoo)
          foo <- foo / sumfoo
          A[, 1, , h] <- lscale + log(foo) # lalpha [ ,1]

          for (T in 2:t) {
            foo <- foo %*% PI[, , h] * f[, T, , h]
            sumfoo <- rowSums(foo)
            lscale <- lscale + log(sumfoo)
            foo <- foo / sumfoo
            A[, T, , h] <- log(foo) + lscale
          }

          llk[h] <- sum(lscale) # log-likelihood

          B[, t, , h] <- 0

          foo <- matrix(rep(rep(1 / k, k), each = num0), ncol = k)
          lscale <- log(k)

          for (T in (t - 1):1) {
            foo <- t(PI[, , h] %*% t(foo * f[, T + 1, , h]))
            B[, T, , h] <- log(foo) + lscale
            sumfoo <- rowSums(foo)
            foo <- foo / sumfoo
            lscale <- lscale + log(sumfoo)
          }

          if (any(prior[h, ] <= 0.05)) {
            evs[h] <- 1
          }
          if (any(is.infinite(A[, , , h]))) {
            llk[h] <- NA
          }
          if (any(is.infinite(B[, , , h]))) {
            llk[h] <- NA
          }
        },
        error = function(e) {
          skip_to_next <<- TRUE
        }
      )

      if (skip_to_next) {
        next
      }
    }

    if (all(evs == 1) == FALSE) {
      llk[evs == 1] <- NA
    }

    df <- data.frame(llk = llk, pos = c(1:nstartR))
    df <- tidyr::drop_na(df)
    df <- df[!is.infinite(rowSums(df)), ]
    bestR <- utils::head(data.table::setorderv(df, cols = "llk", order = -1), n = 1)$pos

    M <- array(M[, , , bestR], dim = c(p, r, k))
    WR <- array(sigmaUY[, , , bestR], dim = c(p, p, k))
    WC <- array(sigmaVY[, , , bestR], dim = c(r, r, k))
    PI2 <- matrix(PI[, , bestR], k, k)
    ef <- array(f[, , , bestR], c(num0, t, k))
    a <- array(A[, , , bestR], c(num0, t, k))
    b <- array(B[, , , bestR], c(num0, t, k))
    lik <- llk[bestR]

    return(list(
      dens = density, model = c(mod.row, mod.col),
      prior = prior[bestR, ], M = M, sigmaUY = WR, sigmaVY = WC, nu = nu, alpha = alpha, eta = eta,
      PI = PI2, f = ef, A = a, B = b, llk = lik
    ))
  }

  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }

  if (any(mod.row == "all")) {
    mod.row.i <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV")
  } else {
    mod.row.i <- mod.row
  }
  if (any(mod.col == "all")) {
    mod.col.i <- c("II", "EI", "VI", "EE", "VE", "EV", "VV")
  } else {
    mod.col.i <- mod.col
  }

  req.model <- expand.grid(mod.row.i, mod.col.i)
  names(req.model) <- paste(rep("V", each = ncol(req.model)), 1:ncol(req.model), sep = "")

  nest.EII.Y <- c("EII", "VII")
  nest.EEI.Y <- c("EEI", "VEI", "EVI", "VVI")
  nest.EEE.Y <- c("EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")
  nest.II.Y <- c("II")
  nest.EI.Y <- c("EI", "VI")
  nest.EE.Y <- c("EE", "VE", "EV", "VV")

  list.comb <- data.frame(matrix(NA, nrow = nrow(req.model), ncol = 2))
  for (i in 1:nrow(req.model)) {
    if (req.model[i, 1] %in% nest.EII.Y) {
      list.comb[i, 1] <- "EII"
    }
    if (req.model[i, 1] %in% nest.EEI.Y) {
      list.comb[i, 1] <- "EEI"
    }
    if (req.model[i, 1] %in% nest.EEE.Y) {
      list.comb[i, 1] <- "EEE"
    }

    if (req.model[i, 2] %in% nest.II.Y) {
      list.comb[i, 2] <- "II"
    }
    if (req.model[i, 2] %in% nest.EI.Y) {
      list.comb[i, 2] <- "EI"
    }
    if (req.model[i, 2] %in% nest.EE.Y) {
      list.comb[i, 2] <- "EE"
    }
  }

  list.comb2 <- unique(list.comb)

  oper <- vector(mode = "list", length = length(k))

  for (g in 1:length(k)) {
    if (verbose == TRUE) {
      print(paste(paste("Initializing Parsimonious", density), "HMMs with k =", k[g]))
    }

    cluster <- snow::makeCluster(nThreads, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)

    pb <- progress::progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = nrow(list.comb2),
      complete = "=", # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">", # Current bar character
      width = 100
    )
    progress <- function(n) {
      pb$tick()
    }
    opts <- list(progress = progress)

    l <- 0

    oper[[g]] <- foreach::foreach(l = 1:nrow(list.comb2), .combine = "comb", .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- r_Pars_init(
        Y = Y, k = k[g], density = density, mod.row = list.comb2[l, 1], mod.col = list.comb2[l, 2], nstartR = nstartR, seed = seed
      )

      list(res)
    }

    if (g == 1) {
      oper2 <- foreach::foreach(i = 1:nrow(req.model), .combine = "comb", .multicombine = TRUE, .init = list(list())) %dopar% {
        for (j in 1:nrow(list.comb2)) {
          if (all(list.comb[i, ] == oper[[1]][[1]][[j]][["model"]])) {
            res <- j
          }
        }

        list(res)
      }
    }

    snow::stopCluster(cluster)
    foreach::registerDoSEQ()
  }

  return(list(
    results = oper,
    k = k,
    req.model = req.model,
    init.used = list.comb,
    index = unlist(oper2[[1]]),
    dens = density
  ))
} # initializing function
