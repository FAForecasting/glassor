source("main.R")

glassof <- function (s, rho, threshold=1e-4, maxit=1e4, approx=FALSE, penalize.diagonal=TRUE, trace=FALSE) {
  n <- nrow(s)
  warm.start <- 0
  ww <- matrix(0,nrow=n,ncol=n)
  xx <- matrix(0,nrow=n,ncol=n)

  rho <- matrix(rho, n,n)

  res <- glassor(n, s, rho, approx, warm.start, trace, penalize.diagonal, threshold, maxit, ww, xx)

  return(res)

}

set.seed(124)
mat <- matrix(rnorm(30*14), 14)
varm <- var(mat)
rho <- 0.1

glassor_res <- glassof(varm, rho)
glasso_res <- glasso::glasso(varm, rho)
huge_res <- huge::huge(varm, lambda = rho, method = "glasso", cov.output = T, verbose = F)

c(max(glassor_res$www - glasso_res$w), max(glasso_res$w - huge_res$cov[[1]]))
c(max(glassor_res$wwwi - glasso_res$wi), max(glasso_res$wi - huge_res$icov[[1]]))



xx <- matrix(c(1,0.4,0.4,0.6), 2, 2)

inv(2, xx, matrix(c(0.3,0.2), 1, 2), xx)

?glasso



options(warn=2)