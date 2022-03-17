source("main.R")

glassof <- function (s, rho, thr=1.0e-4, maxit=1e4, approx=FALSE, penalize.diagonal=TRUE, trace=FALSE) {
  n <- nrow(s)
  warm.start <- 0
  ww <- matrix(0,nrow=n,ncol=n)
  xx <- matrix(0,nrow=n,ncol=n)

  rho <- matrix(rho, n,n)

  res <- glassor(n, s, rho, approx, warm.start, trace, penalize.diagonal, thr, maxit, ww, xx, 1, 1)

  return(res)

}

set.seed(124)
mat <- matrix(rnorm(3*14), 14)

varm <- var(mat)

res2 <- glasso::glasso(varm, 0.1)
res <- glassof(varm, 0.1)

res$www - res2$w
res$wwwi - res2$wi

solve(res2$w) - res2$wi

xx <- matrix(c(1,0.4,0.4,0.6), 2, 2)

inv(2, xx, matrix(c(0.3,0.2), 1, 2), xx)

?glasso



options(warn=2)