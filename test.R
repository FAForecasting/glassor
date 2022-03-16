n <- 3
s <- matrix(c(0.68666584, 0.11700394,  0.07177829, 0.11700394,  1.07891450, -0.09281059, 0.07177829, -0.09281059,  1.04010751), 3, 3, T)
rho <- matrix(0.5, 3,3)
ia <- 0
is <- 0
itrace <- 0
ipen <- 1
thr <- 1e-04
maxit <- 10000
ww <- matrix(0, 3,3)
xx <- ww
niter <- 1
del <- 1
ierr <- 1

glassor(n,s,rho,ia,is,itrace,ipen,thr,maxit,ww,xx,niter,del)

library(glasso)

junk<-.Fortran("glasso",
as.integer(n),
s,
rho,
as.integer(ia),
as.integer(is),
as.integer(itrace),
as.integer(ipen),
thr,
maxit=as.integer(maxit),
ww=ww,
xx=xx,
niter=integer(1),
del=double(1),
ierr=integer(1),
PACKAGE="glasso"
)

