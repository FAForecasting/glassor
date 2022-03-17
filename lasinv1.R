lasinv1 <- function (n, S, rho, approx, warm.start, trace, penalize.diag, threshold, maxit, W, W.inverse, iter) {
  S <- matrix(S, n, n)
  W <- matrix(W, n, n)
  W.inverse <- matrix(W.inverse, n, n)
  rho <- matrix(rho, n, n)

  eps <- 1e-7
  nm1 <- n - 1
  vv <- matrix(0, nm1, nm1)
  if (!approx) xs <- matrix(0, nm1, n)
  s <- numeric(nm1)
  x <- numeric(nm1)
  mm <- numeric(nm1)
  ro <- numeric(nm1)

  #10291
  # Calculate the sum of the absolute values for S^-diag
  shr <- 0
  for (j in 1:n) {
    for (k in 1:n) {
      if (j != k) {
        shr <- shr + abs(S[j, k])
      }
    }
  }

  #10301
  if (shr == 0) {
    W[,] <- 0
    W.inverse[,] <- 0
      for (j in 1:n) {
        if (penalize.diag == 0) {
          W[j, j] <- S[j, j]
        } else {
          W[j, j] <- S[j, j] + rho[j, j]
        }
        W.inverse[j, j] <- 1.0 / max(W[j, j], eps)
      }
      return(list(ww=W, wwi=W.inverse, del=0, iter=iter))
  }
  #10331
  shr <- threshold * shr / nm1

  if (approx) {
    if (warm.start == 0) W.inverse[,] <- 0
    for (m in 1:n) {
      res <- setup(m, n, S, rho, S, vv, s, ro)
      vv <- res$vv
      s <- res$s
      ro <- res$r
      l <- 0
      for (j in 1:n) {
        if (j != m) {
          l <- l+1
          x[l] <- W.inverse[j, m]
        }
      }
      res <- lasso(ro, nm1, vv, s, shr/n, x, mm)
      s <- res$s
      x <- res$x
      mm <- res$mm
      l <- 0
      for (j in 1:n) {
        if (j != m) {
          l <- l+1
          W.inverse[j, m] <- x[l]
        }
      }
    }
    #10401
    iter <- 1
    return(list(ww=W, wwi=W.inverse, del=del, iter=iter))
  }

  if (!warm.start) {
    W <- S
    xs[,] <- 0
  } else {
    for (j in 1:n) {
      xjj <- -W.inverse[j, j]
      l <- 0
      for (k in 1:n) {
        if (k != j) {
          l <- l + 1
          xs[l,j] <- W.inverse[k, j] / xjj
        }
      }
    }
  }
  #10451

  for (j in 1:n) {
    if (penalize.diag == 0) {
      W[j, j] <- S[j, j]
    } else {
      W[j, j] <- S[j, j] + rho[j, j]
    }
  }

  #10481
  iter <- 0
  W.delta <- Inf
  while (W.delta >= shr && iter < maxit) {
    for (m in 1:n) {
      if (trace) {
        cat("m: ", m, "\n")
      }
      x <- xs[, m]
      ws <- W[, m]
      res <- setup(m, n, S, rho, W, vv, s, ro)
      vv <- res$vv
      s <- res$s
      ro <- res$r
      so <- s
      res <- lasso(ro, nm1, vv, s, shr / sum(abs(vv)), x, mm)
      s <- res$s
      x <- res$x
      mm <- res$mm
      l <- 0
      for (j in 1:n) {
        if (j != m) {
          l <- l + 1
          W[j, m] <- so[l] - s[l]
          W[m, j] <- W[j, m]
        }
      }
      W.delta <- sum(abs(W[, m]-ws))
      xs[,m] <- x
    }

    iter <- iter + 1
  }
  del <- W.delta / nm1
  W.inverse <- inv(n, W, xs, W.inverse)
  return(list(ww=W, wwi=W.inverse, del=del, iter=iter))
}
