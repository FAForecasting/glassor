#' Select the best result from a path estimate based on ebic
#' @param results The results from glassorpath
#' @return The result with the lowest ebic
#' @export
glassor.select <- function (results) {
  ebic <- lapply(results, ebic.glassor)
  best <- which.min(ebic)
  return(results[[best]])
}

#' Calculate the ebic for a result
#' @param res The result
#' @return The ebic of res
#' @export
ebic.glassor <- function (res) {
  df <- sum(res$wi[lower.tri(res$wi)] != 0)
  d <- ncol(res$w)
  ebic <- -2 * res$loglik + log(res$nobs) * df + 4 * res$rho[1] * log(d) * df
  return(ebic)
}
