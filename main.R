

glassor <- function (nn, s, rrho, approx, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv) {
  if (approx) {
    res <- lasinv1(nn, s, rrho, approx, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv, 1)
    return(res)
  } else {
    # 10021
    res <- connect(nn, s, rrho)
    nc <- res$nc
    ic <- res$ic
    ir <- res$ir

    nnq <- 0
    for (kc in 1:nc) {
      nnq <- max(ic[2,kc] - ic[1,kc]+1, nnq)
    }
    # 10031
    nnq <- nnq^2
    ss <- numeric(nnq)
    rho <- numeric(nnq)
    ww <- numeric(nnq)
    wwi <- numeric(nnq)

    niter <- 0
    ddel <- 0
    l <- 0
    #10041
    for (kc in 1:nc) {
      n <- ic[2, kc] - ic[1, kc]+1
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
        res <- lasinv1(n, ss, rho, approx, warm.start, trace, penalize.diag, threshold, maxit, ww, wwi, niter)
        ww <- res$ww
        wwi <- res$wwi
        niter <- res$niter + 1
        ddel <- ddel + res$del
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
    ddel <- ddel/nc
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
  }
  return(list(www=W, wwwi=W.inv, del=ddel, niter=niter))
}



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
        iss <- is
        il <- is + res$na - 1
        if (il >= n) break
        is <- is + res$na
        res <- row(nc, res$na, ir[iss:n], n, S, rho, ie, ir[is:n])
        ie <- res$ie
        ir[is:n] <- res$ir
      }
      #10252
      ic[2, nc] <- il
    }
  }
  return(list(nc=nc, ic=ic, ir=ir))
}

row <- function (nc, nr, jr, n, ss, rho, ie, kr) {
  na <- 0
  for (l in 1:nr) {
    k <- jr[l]
    for (j in 1:n) {
      if (ie[j] <= 0 && j != k && abs(ss[j, k]) > rho[j, k]) {
        na <- na + 1
        kr[na] <- j
        ie[j] <- nc
      }
    }
  }
  return(list(ie=ie, na=na, ir=kr))
}

lasinv1 <- function (n, S, rho, approx, warm.start, trace, penalize.diag, threshold, maxit, ww, wwi, niter) {
  S <- matrix(S, n, n)
  ww <- matrix(ww, n, n)
  wwi <- matrix(wwi, n, n)
  rho <- matrix(rho, n, n)

  eps <- 1e-7
  nm1 <- n - 1
  vv <- matrix(0, nm1, nm1)
  if (!approx) xs <- matrix(0, nm1, n)
  s <- numeric(nm1)
  x <- numeric(nm1)
  z <- numeric(nm1)
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
    ww[,] <- 0
    wwi[,] <- 0
      for (j in 1:n) {
        if (penalize.diag == 0) {
          ww[j,j] <- S[j, j]
        } else {
          ww[j,j] <- S[j, j] + rho[j, j]
        }
        wwi[j,j] <- 1.0 / max(ww[j,j], eps)
      }
      return(list(ww=ww, wwi=wwi, del=0, niter=niter))
  }
  #10331
  shr <- threshold * shr / nm1

  if (approx) {
    if (warm.start == 0) wwi[,] <- 0
    for (m in 1:n) {
      res <- setup(m, n, S, rho, S, vv, s, ro)
      vv <- res$vv
      s <- res$s
      ro <- res$r
      l <- 0
      for (j in 1:n) {
        if (j != m) {
          l <- l+1
          x[l] <- wwi[j,m]
        }
      }
      res <- lasso(ro, nm1, vv, s, shr/n, x, z, mm)
      s <- res$s
      x <- res$x
      z <- res$z
      mm <- res$mm
      l <- 0
      for (j in 1:n) {
        if (j != m) {
          l <- l+1
          wwi[j,m] <- x[l]
        }
      }
    }
    #10401
    niter <- 1
    return(list(ww=ww, wwi=wwi, del=del, niter=niter))
  }
  if (!warm.start) {
    ww <- S
    xs[,] <- 0
  } else {
    for (j in 1:n) {
      xjj <- -wwi[j,j]
      l <- 0
      for (k in 1:n) {
        if (k != j) {
          l <- l + 1
          xs[l,j] <- wwi[k, j] / xjj
        }
      }
    }
  }
  #10451

  for (j in 1:n) {
    if (penalize.diag == 0) {
      ww[j,j] <- S[j, j]
    } else {
      ww[j,j] <- S[j, j] + rho[j, j]
    }
  }

  #10481
  niter <- 0
  dlx <- 0
  while (dlx < shr && niter < maxit) {
    dlx <- 0
    for (m in 1:n) {
      if (trace) {
        cat("m: ", m, "\n")
      }
      x <- xs[, m]
      ws <- ww[, m]
      res <- setup(m, n, S, rho, ww, vv, s, ro)
      vv <- res$vv
      s <- res$s
      ro <- res$r
      so <- s
      res <- lasso(ro, nm1, vv, s, shr / sum(abs(vv)), x, z, mm)
      s <- res$s
      x <- res$x
      z <- res$z
      mm <- res$mm
      l <- 0
      for (j in 1:n) {
        if (j != m) {
          l <- l + 1
          ww[j,m] <- so[l] - s[l]
          ww[m,j] <- ww[j,m]
        }
      }
      dlx <- max(dlx, sum(abs(ww[,m]-ws)))
      xs[,m] <- x
    }

    niter <- niter + 1
  }
  del <- dlx / nm1
  wwi <- inv(n, ww, xs, wwi)
  return(list(ww=ww, wwi=wwi, del=del, niter=niter))
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

lasso <- function (rho, n, vv, s, threshold, x, z, mm) {
  res <- fatmul(2, n, vv, x, s, z, mm)
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
  return(list(s=s, x=x, z=z, mm=mm))
}

fsign <- function (x, y) {
  x <- abs(x)
  if (y > 0) x
  else -x
}

fatmul <- function (it, n, vv, x, s, z, m) {
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

