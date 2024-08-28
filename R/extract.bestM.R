#' Selection of the best fitting model(s)
#'
#' This functions extracts the best fitting model(s) according to the Bayesian information criterion (BIC).
#'
#' @param results The output of the \code{Eigen.HMM_fit()} function.
#' @param top Integer. Specifies the number of top-ranked models to display based on the Bayesian Information Criterion (BIC).
#'
#' @return A list containing the required best fitting model(s).
#' @export
#' @importFrom foreach %dopar%
#' @examples
#' data(simData)
#' Y <- simData$Y
#' init <- Eigen.HMM_init(Y = Y, k = 2, density = "MVT", mod.row = "EEE", mod.col = "EE", nstartR = 10)
#' fit <- Eigen.HMM_fit(Y = Y, init.par = init, nThreads = 1)
#' win <- extract.bestM(results = fit, top = 1)
extract.bestM <- function(results, top = 1) {
  k <- length(results[["results"]])
  num.mod <- length(results[["results"]][[1]][[1]])
  list.mod <- results[["models"]]
  list.mod2 <- do.call("rbind", replicate(k, list.mod, simplify = FALSE))
  count.k <- sort(rep(1:k, num.mod))
  count.mod <- rep(1:num.mod, k)
  list.mod3 <- data.frame(list.mod2, count.k, count.mod)

  allBIC <- numeric(k * num.mod)

  cont <- sp.th <- 0

  for (j in 1:k) {
    for (i in 1:num.mod) {
      if (!all(is.na(results[["results"]][[j]][[1]][[i]]))) {
        cont <- cont + 1
        allBIC[cont] <- -results[["results"]][[j]][[1]][[i]][["BIC"]]

        tclass <- results[["results"]][[j]][[1]][[i]][["class"]]
        t <- ncol(tclass)
        kk <- length(results[["results"]][[j]][[1]][[i]][["prior"]])
        sp.mat <- matrix(0, nrow = kk, ncol = t)
        for (l in 1:t) {
          if (any(table(tclass[, l]) / sum(table(tclass[, l])) <= sp.th) == TRUE) {
            sp.mat[as.numeric(names(which(table(tclass[, l]) / sum(table(tclass[, l])) <= sp.th))), l] <- 1
          }
        }

        f.ev <- rowSums(sp.mat)
        if (any(f.ev == t)) {
          allBIC[cont] <- NA
        }
      } else {
        cont <- cont + 1
        allBIC[cont] <- NA
      }
    }
  }

  topBIC <- which(allBIC >= sort(allBIC, decreasing = T)[top])
  topBIC.order <- order(allBIC[topBIC], decreasing = T)
  tempBIC <- list.mod3[topBIC[topBIC.order], ]
  bestBIC <- vector(mode = "list", length = top)
  for (i in 1:top) {
    bestBIC[[i]] <- results[["results"]][[as.numeric(tempBIC[i, 3])]][[1]][[as.numeric(tempBIC[i, 4])]]
  }

  return(bestBIC = bestBIC)
} # extract best fitting model
