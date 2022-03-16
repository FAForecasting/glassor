glassof <- function (s, rho, thr=1.0e-4, maxit=1e4, approx=FALSE, penalize.diagonal=TRUE, trace=FALSE) {
  n <- nrow(s)
  warm.start <- 0
  ww <- matrix(0,nrow=n,ncol=n)
  xx <- matrix(0,nrow=n,ncol=n)

  rho <- matrix(rho, n,n)

  res <- glassor(n, s, rho, approx, warm.start, itrace, ipen, thr, maxit, ww, xx, 1, 1)

  return(res)

}

glassof()

mat <- matrix(rnorm(3*14), 14)

varm <- var(mat)

res2 <- glasso::glasso(varm, 0.2)
res <- glassof(varm, 0.2)

max(res$www - res2$w)
res$wwwi - res2$wi

solve(res2$w) - res2$wi

inv(3, res2$w, matrix(0, 3, 3))

?glasso



options(warn=2)