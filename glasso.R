source("main.R")



set.seed(124)
mat <- matrix(rnorm(30*14), 14)
varm <- var(mat)
rho <- 0.1

glassor_res <- glassof(varm, rho)
glasso_res <- glasso::glasso(varm, n=14, rho)
huge_res <- huge::huge(varm, lambda = rho, method = "glasso", cov.output = T, verbose = F)

c(max(glassor_res$www - glasso_res$w), max(glasso_res$w - huge_res$cov[[1]]))
c(max(glassor_res$wwwi - glasso_res$wi), max(glasso_res$wi - huge_res$icov[[1]]))



xx <- matrix(c(1,0.4,0.4,0.6), 2, 2)

inv(2, xx, matrix(c(0.3,0.2), 1, 2), xx)

?glasso



options(warn=2)