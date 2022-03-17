#' Select the best result from a path estimate based on ebic
#' @param results The results from glassorpath
#' @return The result with the lowest ebic
#' @export
glassor.select <- function (results) {
  ebic <- lapply(results, ebic.glasso)
  best <- which.min(ebic)
  return(results[[best]])
}

#' Calculate the ebic for a result
#' @param res The result
#' @return The ebic of res
#' @export
ebic.glasso <- function (res) {
  d <- ncol(res$w)
  ebic <- -res$nobs * res$loglik + log(res$nobs) * res$df + 4 * ebic.gamma * log(d) * res$df
  return(ebic)
}