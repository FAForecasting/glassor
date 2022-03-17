glassof <- function (s, rho, threshold=1e-4, maxit=1e4, approx=FALSE, penalize.diagonal=TRUE, trace=FALSE) {
  n <- nrow(s)
  warm.start <- 0
  ww <- matrix(0,nrow=n,ncol=n)
  xx <- matrix(0,nrow=n,ncol=n)

  rho <- matrix(rho, n,n)

  if (approx) res <- glassor.approx(n, s, rho, warm.start, trace, penalize.diagonal, threshold, maxit, ww, xx)
  else res <- glassor(n, s, rho, warm.start, trace, penalize.diagonal, threshold, maxit, ww, xx)


  return(res)

}

glassor <- function (nn, s, rrho, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv) {
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

  iter <- 0
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
      res <- lasinv1(n, ss, rho, FALSE, warm.start, trace, penalize.diag, threshold, maxit, ww, wwi, iter)
      ww <- res$ww
      wwi <- res$wwi
      iter <- res$iter + 1
      del <- res$del
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

  return(list(www=W, wwwi=W.inv, del=res$del/nc, iter=iter))
}

glassor.approx <- function (nn, s, rrho, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv) {
  lasinv1(nn, s, rrho, TRUE, warm.start, trace, penalize.diag, threshold, maxit, W, W.inv, 1)
}
