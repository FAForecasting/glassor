#' Graphical lasso, based on the Fortran version by Tibshirani et. al.
#' @param s The observed covariance matrix
#' @param rho The regularisation parameter for the lasso
#' @param nobs The number of observations, needed to calculate the log likelihood
#' @export
glassor <- function (s, rho, nobs=NULL, threshold=1e-4, maxit=1e4,
                     approx=FALSE, penalize.diagonal=TRUE,
                     start=c("cold","warm"), w.init=NULL, wi.init=NULL, trace=FALSE) {

  # Set up the algorithm for cold or warm start
  n <- nrow(s)
  start.type <- match.arg(start)
  if (start.type == "cold") {
    warm.start <- 0
    W <- matrix(0,nrow=n,ncol=n)
    W.inv <- matrix(0,nrow=n,ncol=n)
  } else if (start.type == "warm") {
    warm.start <- 1
    if (is.null(w.init) | is.null(wi.init)) {
      stop("Warm start specified: w.init and wi.init must be non-null")
    }
    W <- w.init
    W.inv <- wi.init
  }

  # Make sure that rho has the right form, it can be supplied both as a double and a matrix.
  rho_mat <- matrix(rho, n,n)

  # Run the algorithm, exact or approximate
  if (approx) res <- .glasso.approx(n, s, rho_mat, warm.start, trace, penalize.diagonal, threshold, maxit, W, W.inv)
  else res <- .glasso(n, s, rho_mat, warm.start, trace, penalize.diagonal, threshold, maxit, W, W.inv)

  # Calculate log likelihood if nobs is set
  if (!is.null(nobs)) res$loglik <- loglik(s, rho_mat, nobs, res$wi, penalize.diagonal)
  else res$loglik <- NA
  res$nobs <- nobs
  res$rho <- rho

  return(res)
}

#' Graphical lasso, based on the Fortran version by Tibshirani et. al.
#' @param s The observed covariance matrix
#' @param rho The regularisation parameter for the lasso
#' @param nobs The number of observations, needed to calculate the log likelihood
#' @export
glassorpath <- function (s, rho, nobs=NULL, threshold=1e-4, maxit=1e4,
                         approx=FALSE, penalize.diagonal=TRUE,
                         start=c("cold","warm"), w.init=NULL, wi.init=NULL, trace=FALSE, parallel=FALSE) {
  if (is.vector(rho)) {
    rho <- as.list(rho)
  } else if (is.matrix(rho)) {
    rho <- list(rho)
  }

  if (!parallel) runfun <- lapply
  else runfun <- mclapply

  results <- runfun(rho, glassor, s=s, nobs=nobs, threshold=threshold, maxit=maxit,
                                   approx=approx, penalize.diagonal=penalize.diagonal,
                                   start=start, w.init=w.init, wi.init=wi.init, trace=trace)
  return(results)
}

loglik <- function (s, rho, nobs, W.inverse, penalize.diagonal) {
  W.inv.det <- det(W.inverse)

  crit <- -log(W.inv.det) + sum(diag(s %*% W.inverse))
  if (!penalize.diagonal) {
    diag(W.inverse) <- 0
  }
  crit <- crit + sum(abs(rho * W.inverse))

  loglik <- -(nobs / 2) * crit

  return(loglik)
}

.glasso <- function (nn, s, rrho, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv) {
  # 10021
  res <- connect(nn, s, rrho)
  nc <- res$nc
  ic <- res$ic
  ir <- res$ir

  nnq <- 0
  for (kc in 1:nc) {
    nnq <- max(ic[2,kc] - ic[1,kc] + 1, nnq)
  }
  # 10031
  nnq <- nnq^2
  ss <- numeric(nnq)
  rho <- numeric(nnq)
  ww <- numeric(nnq)
  wwi <- numeric(nnq)

  iter <- 0
  l <- 0
  #10041
  for (kc in 1:nc) {
    n <- ic[2, kc] - ic[1, kc] + 1
    if (n <= 1) {
      k <- ir[ic[1, kc]]
      W[, k] <- 0
      W[k,] <- 0
      W.inv[, k] <- 0
      W.inv[k,] <- 0
    } else {
      # 10061
      kbegin <- ic[1, kc]
      kend <- ic[2, kc]
      l <- 0
      for (k in kbegin:kend) {
        ik <- ir[k]
        for (j in kbegin:kend) {
          ij <- ir[j]
          l <- l+1
          ss[l] <- s[ij, ik]
          rho[l] <- rrho[ij,ik]
          ww[l] <- W[ij, ik]
          wwi[l] <- W.inv[ij, ik]
        }
        #10081
      }
      #10071
      res <- lasinv1(n, ss, rho, FALSE, warm.start, trace, penalize.diag, threshold, maxit, ww, wwi, iter)
      ww <- res$ww
      wwi <- res$wwi
      iter <- res$iter + 1
      del <- res$del
      for (j in kbegin:kend) {
        k <- ir[j]
        W[, k] <- 0
        W[k,] <- 0
        W.inv[, k] <- 0
        W.inv[k,] <- 0
      }
      #10091
      l <- 0
      for (k in kbegin:kend) {
        ik <- ir[k]
        for (j in kbegin:kend) {
          l <- l+1
          W.inv[ir[j], ik] <- wwi[l]
        }
      }
      #10111
      l <- 0
      for (k in kbegin:kend) {
        ik <- ir[k]
        for (j in kbegin:kend) {
          l <- l+1
          W[ir[j], ik] <- ww[l]
        }
      }
    }
  }
  # 10041
  for (j in 1:nn) {
    if (W[j, j] == 0) {
      if (penalize.diag == 0) {
        W[j, j] <- s[j, j]
      } else {
        W[j, j] <- s[j, j] + rrho[j, j]
      }
      W.inv[j, j] <- 1 / W[j, j]
    }
  }

  return(list(w=W, wi=W.inv, del=res$del/nc, iter=iter))
}

.glasso.approx <- function (nn, s, rrho, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv) {
  lasinv1(nn, s, rrho, TRUE, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv, 1)
}
