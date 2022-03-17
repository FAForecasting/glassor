
connect <- function (n, S, rho) {
  ic <- matrix(0,2, n)
  ir <- numeric(n)
  ie <- numeric(n)
  nc <- 0
  is <- 1
  for (k in 1:n) {
    if (ie[k] <= 0) {
      ir[is] <- k
      nc <- nc + 1
      ie[k] <- nc
      ic[1, nc] <- is
      is <- is + 1

      res <- row(nc, 1, ir[(is-1):n], n, S, rho, ie, ir[is:n])
      ie <- res$ie
      ir[is:n] <- res$ir

      il <- is - 1
      while (res$na != 0) {
        il <- is + res$na - 1
        if (il >= n) break
        is <- is + res$na
        res <- row(nc, res$na, ir[(is-res$na):n], n, S, rho, ie, ir[is:n])
        ie <- res$ie
        ir[is:n] <- res$ir
      }
      #10252
      ic[2, nc] <- il
    }
  }
  return(list(nc=nc, ic=ic, ir=ir))
}

row <- function (nc, nr, jr, n, ss, rho, ie, ir) {
  na <- 0
  for (l in 1:nr) {
    k <- jr[l]
    for (j in 1:n) {
      if (ie[j] <= 0 && j != k && abs(ss[j, k]) > rho[j, k]) {
        na <- na + 1
        ir[na] <- j
        ie[j] <- nc
      }
    }
  }
  return(list(ie=ie, na=na, ir=ir))
}

