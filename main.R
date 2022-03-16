

glassor <- function (nn, s, rrho, approx, warm.start, trace, penalize.diag, threshold, maxit, www, wwwi, niter, ddel) {
  if (approx != 0) {
    res <- lasinv1(nn, s, rrho, approx, warm.start, trace, penalize.diag, threshold, maxit, www, wwwi, niter)
    return(res)
  } else {
    # 10021
    ic <- matrix(0,2, nn)
    ir <- numeric(nn)
    ie <- numeric(nn)
    res <- connect(nn, s, rrho, nc, ic, ir, ie)
    nc <- res$nc
    ic <- res$ic
    ir <- res$ir
    ie <- res$ie
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
        www[,k] <- 0
        www[k,] <- 0
        wwwi[,k] <- 0
        wwwi[k,] <- 0
      } else {
        # 10061
        kb <- ic[1, kc]
        ke <- ic[2, kc]
        l <- 0
        for (k in kb:ke) {
          ik <- ir[k]
          for (j in kb:ke) {
            ij <- ir[j]
            l <- l+1
            ss[l] <- s[ij, ik]
            rho[l] <- rrho[ij,ik]
            ww[l] <- www[ij,ik]
            wwi[l] <- wwwi[ij,ik]
          }
          #10081
        }
        #10071
        res <- lasinv1(n, ss, rho, approx, warm.start, trace, penalize.diag, threshold, maxit, ww, wwi, niter)
        ww <- res$ww
        wwi <- res$wwi
        niter <- res$niter + 1
        ddel <- ddel + res$del
        for (j in kb:ke) {
          k <- ir[j]
          www[,k] <- 0
          www[k,] <- 0
          wwwi[,k] <- 0
          wwwi[k,] <- 0
        }
        #10091
        l <- 0
        for (k in kb:ke) {
          for (j in kb:ke) {
            l <- l+1
            wwwi[ir[j],ik] <- wwi[l]
          }
        }
        #10111
        if (approx == 0) {
          l <- 0
          for (k in kb:ke) {
            ik <- ir[k]
            for (j in kb:ke) {
              l <- l+1
              www[ir[j],ik] <- ww[l]
            }
          }
        }
      }
    }
    # 10041
    ddel <- ddel/nc
    if (approx == 0) {
      for (j in 1:nn) {
        if (www[j,j] == 0) {
          if (penalize.diag == 0) {
            www[j,j] <- s[j, j]
          } else {
            www[j,j] <- s[j, j] + rrho[j, j]
          }
          wwwi[j,j] <- 1 / www[j,j]
        }
      }
    }
  }
  return(list(www=www, wwwi=wwwi, del=ddel, niter=niter))
}



connect <- function (n,ss,rho,nc,ic,ir,ie) {
  ie[seq_along(ie)] <- 0
  nc <- 0
  is <- 1
  for (k in 1:n) {
    if (ie[k] <= 0) {
      ir[is] <- k
      nc <- nc + 1
      ie[k] <- nc
      ic[1,nc] <- is
      is <- is+1
      res <- row(nc,1,ir[(is-1):n],n,ss,rho,ie,na,ir[is:n])
      ie <- res$ie
      na <- res$na
      ir[is:n] <- res$ir
      if (na == 0) {
        ic[2,nc] <- is-1
      } else {
        while (TRUE) {
          nas <- na
          iss <- is
          il <- iss + nas - 1
          if (il > n) break
          is <- is + na
          res <- row(nc,nas,ir[iss:n],n,ss,rho,ie,na,ir[is:n])
          ie <- res$ie
          na <- res$na
          ir[is:n] <- res$ir
          if (na == 0) break
        }
        #10252
        ic[2, nc] <- il
      }
    }
  }

  return(list(nc=nc, ic=ic, ir=ir, ie=ie))
}

row <- function (nc, nr, jr, n, ss, rho, ie, na, kr) {
  na <- 0
  for (l in 1:nr) {
    k <- jr[l]
    for (j in 1:n) {
      if (ie[j] <= 0 && j != k && abs(ss[j,k]) > rho[j,k]) {
        na <- na+1
        kr[na] <- j
        ie[j] <- nc
      }
    }
  }
  return(list(ie=ie, na=na, ir=kr))
}

lasinv1 <- function (n, ss, rho, approx, is, trace, penalize.diag, threshold, maxit, ww, wwi, niter) {
  ss <- matrix(ss, n, n)
  ww <- matrix(ww, n, n)
  wwi <- matrix(wwi, n, n)
  rho <- matrix(rho, n, n)

  eps <- 1e-7
  nm1 <- n - 1
  vv <- matrix(0, nm1, nm1)
  if (approx == 0) xs <- matrix(0, nm1, n)
  s <- numeric(nm1)
  s0 <- numeric(nm1)
  x <- numeric(nm1)
  z <- numeric(nm1)
  mm <- numeric(nm1)
  ro <- numeric(nm1)
  if (approx == 0) {
    ws <- numeric(n)
  }

  #10291
  shr <- 0
  for (j in 1:n) {
    for (k in 1:n) {
      if (j != k) {
        shr <- shr + abs(ss[j,k])
      }
    }
  }

  #10301
  if (shr == 0) {
    ww[,] <- 0
    wwi[,] <- 0
      for (j in 1:n) {
        if (penalize.diag == 0) {
          ww[j,j] <- ss[j,j]
        } else {
          ww[j,j] <- ss[j,j] + rho[j,j]
        }
        wwi[j,j] <- 1.0 / max(ww[j,j], eps)
      }
      return(list(ww=ww, wwi=wwi, del=0, niter=niter))
  }
  #10331
  shr <- threshold * shr / nm1

  if (approx != 0) {
    if (is == 0) wwi[,] <- 0
    for (m in 1:n) {
      res <- setup(m,n,ss,rho,ss,vv,s,ro)
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
  if (is == 0) {
    ww <- ss
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
      ww[j,j] <- ss[j,j]
    } else {
      ww[j,j] <- ss[j,j] + rho[j,j]
    }
  }

  #10481
  niter <- 0
  dlx <- 0
  while (dlx < shr && niter < maxit) {
    for (m in 1:n) {
      if (trace != 0) {
        intpr('m',1,m,1)
      }
      x <- xs[, m]
      ws <- ww[, m]
      res <- setup(m, n, ss, rho, ww, vv, s, ro)
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

setup <- function (m,n,ss,rho,ww,vv,s,r) {
  l <- 0
  for (j in 1:n) {
    if (j != m) {
      l <- l + 1
      r[l] <- rho[j,m]
      s[l] <- ss[j,m]
      i <- 0
      for (k in 1:n) {
        if (k != m) {
          i <- i + 1
          vv[i,l] <- ww[k,j]
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

fatmul <- function (it,n,vv,x,s,z,m) {
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

  if (l > fac*n) {
    if (it == 1) {
      s <- matmul[vv,x]
    } else {
      s <- s - matmul[vv,x]
    }
  } else {
    if (it == 1) {
      for (j in 1:n) {
        s[j] <- dot_product(vv[j, m[1:l]], x[1:l])
      }
    } else {
      for (j in 1:n) {
        s[j] <- s[j] - dot_product(vv[m[1:l],j],z[1:l])
      }
    }
  }
  return(list(s=s, z=z, m=m))
}

inv <- function (n, ww, xs, wwi) {
  nm1 <- n - 1
  xs <- -xs
  wwi[1,1] <- 1 / (ww[1,1] + dot_product(xs[,1],ww[2:n,1]))
  wwi[2:n,1] <- wwi[1,1] * xs[,1]
  wwi[n,n]  <- 1/ (ww[n,n] + dot_product(xs[,n],ww[1:nm1,n]))
  wwi[1:nm1,n] <- wwi[n,n] * xs[,n]

  for (j in 2:nm1) {
    jm1 <- j-1
    jp1 <- j+1
    wwi[j,j] <- 1 / (ww[j,j] + dot_product(xs[1:jm1,j], ww[1:jm1,j]) + dot_product(xs[j:nm1,j],ww[jp1:n,j]))
    wwi[1:jm1,j] <- wwi[j,j] * xs[1:jm1,j]
    wwi[jp1:n,j] <- wwi[j,j] * xs[j:nm1,j]
  }
  return(wwi)
}

dot_product <- function (x,y) {
  sum(x * y)
}

