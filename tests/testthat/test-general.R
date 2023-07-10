test_that("Tests to that the results are equivalent to glasso", {
  set.seed(12345)
  for (i in 2:20) {
    for (j in 1:10) {
      nobs <- sample.int(100, 1) + i
      data <- matrix(rnorm(i*nobs), nobs)
      s_observed <- var(data)
      rho <- j / 10

      glassor_res <- glassor(s_observed, rho, nobs=nobs)
      glasso_res <- glasso::glasso(s_observed, rho, nobs=nobs)

      expect_equal(glassor_res$w, glasso_res$w, tolerance = 8e-7)
      expect_equal(glassor_res$wi, glasso_res$wi, tolerance = 8e-7)
      expect_equal(glassor_res$loglik / glasso_res$loglik, 1, tolerance = 5e-7)

      cat(j)

    }
    cat(i)
  }
})

test_that("Tests glasso path, single and multi threaded", {
  set.seed(12345)
  for (i in 5:10) {
    for (j in 1:10) {
      nobs <- sample.int(100, 1) + i*2
      data <- matrix(rnorm(i*2*nobs), nobs)
      s_observed <- var(data)
      rho <- seq(j / 10, 1, 0.1)

      glassor_res <- glassorpath(s_observed, rho, nobs=nobs)
      glassor_res <- glassorpath(s_observed, rho, nobs=nobs, parallel=TRUE)

      select_res <- glassor.select(glassor_res)


      expect_equal(glassor_res$w, glasso_res$w, tolerance = 8e-3)
      expect_equal(glassor_res$wi, glasso_res$wi, tolerance = 8e-3)
      expect_equal(glassor_res$loglik / glasso_res$loglik, 1, tolerance = 5e-3)

      cat(j)

    }
    cat(i)
  }
})
