source("lasinv1.R")


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


setup <- function (m, n, S, rho, W, vv, s, r) {
  l <- 0
  for (j in 1:n) {
    if (j != m) {
      l <- l + 1
      r[l] <- rho[j,m]
      s[l] <- S[j, m]
      i <- 0
      for (k in 1:n) {
        if (k != m) {
          i <- i + 1
          vv[i,l] <- W[k, j]
        }
      }
    }
  }
  return(list(vv=vv, s=s, r=r))
}

lasso <- function (rho, n, vv, s, threshold, x, mm) {
  res <- fatmul(2, n, vv, x, s, mm)
  s <- res$s
  z <- res$z
  mm <- res$mm
  while (TRUE) {
    dlx <- 0
    for (j in 1:n) {
      xj <- x[j]
      x[j] <- 0.0
      t <- s[j] + vv[j,j] * xj
      if (abs(t)-rho[j] > 0) x[j] <- fsign(abs(t) - rho[j],t) / vv[j,j]
      if (x[j] != xj) {
        del <- x[j] - xj
        dlx <- max(dlx, abs(del))
        s <- s - del * vv[,j]
      }
    }
    if (dlx < threshold) break
  }
  return(list(s=s, x=x, mm=mm))
}

fsign <- function (x, y) {
  x <- abs(x)
  if (y > 0) x
  else -x
}

fatmul <- function (it, n, vv, x, s, m) {
  z <- numeric(n)
  fac <- 0.2

  l <- 0
  for (j in 1:n) {
    if (x[j] != 0) {
      l <- l + 1
      m[l] <- j
      z[l] <- x[j]
    }
  }

  #10611

  if (l > fac * n) {
    if (it == 1) {
      s <- vv %*% x
    } else {
      s <- s - vv %*% x
    }
  } else {
    if (it == 1) {
      for (j in 1:n) {
        s[j] <- dot_product(vv[j, m[1:l]], z[1:l])
      }
    } else {
      for (j in 1:n) {
        s[j] <- s[j] - dot_product(vv[m[1:l], j], z[1:l])
      }
    }
  }
  return(list(s=s, z=z, m=m))
}

inv <- function (n, ww, xs, wwi) {
  nm1 <- n - 1
  xs <- -xs
  wwi[1,1] <- 1 / (ww[1,1] + dot_product(xs[,1], ww[2:n,1]))
  wwi[2:n,1] <- wwi[1,1] * xs[,1]
  wwi[n,n]  <- 1/ (ww[n,n] + dot_product(xs[,n], ww[1:nm1,n]))
  wwi[1:nm1,n] <- wwi[n,n] * xs[,n]

  if (nm1 >= 2) {
    for (j in 2:nm1) {
      jm1 <- j-1
      jp1 <- j+1
      wwi[j,j] <- 1 / (ww[j,j] + dot_product(xs[1:jm1,j], ww[1:jm1,j]) + dot_product(xs[j:nm1,j],ww[jp1:n,j]))
      wwi[1:jm1,j] <- wwi[j,j] * xs[1:jm1,j]
      wwi[jp1:n,j] <- wwi[j,j] * xs[j:nm1,j]
    }
  }
  return(wwi)
}

dot_product <- function (x,y) {
  sum(x * y)
}

