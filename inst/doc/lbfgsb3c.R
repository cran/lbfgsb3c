## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
## Added 20190717 to get vignette to build
pkgbuild::compile_dll()

## ---- bt, echo=TRUE-----------------------------------------------------------
# ref BT.RES in Nash and Walker-Smith (1987)

bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

bt.badsetup<-function(n){
   x<-rep(0,n)
   lo<-rep(0,n)
   up<-lo # to get arrays set
   bmsk<-rep(1,n)
   bmsk[(trunc(n/2)+1)]<-0
   for (i in 1:n) { 
      x[i]<-2.2*i-n
      lo[i]<-1.0*(i-1)*(n-1)/n
      up[i]<-1.0*i*(n+1)/n
   }
   result<-list(x=x, lower=lo, upper=up, bdmsk=bmsk)
}

bt.setup0<-function(n){
   x<-rep(0,n)
   lo<-rep(0,n)
   up<-lo # to get arrays set
   bmsk<-rep(1,n)
   bmsk[(trunc(n/2)+1)]<-0
   for (i in 1:n) { 
      lo[i]<-1.0*(i-1)*(n-1)/n
      up[i]<-1.0*i*(n+1)/n
   }
   x<-0.5*(lo+up)
   result<-list(x=x, lower=lo, upper=up, bdmsk=bmsk)
}

## library(lbfgsb3c)
library(optimx)
# sessionInfo()
oppv<- packageVersion("optimx")
cat("optimx version is ", as.character(oppv),"\n")
if (as.character(oppv) < '2019-6.14' ) {
  cat("Note -- Older version of optimx in use\n")
meths <- c("L-BFGS-B", "Rvmmin", "Rcgmin", "Rtnmin")
nn <- 4
goody <- bt.setup0(nn)
lo <- goody$lower
up <- goody$upper
x0 <- goody$x
goody
## optim()
## Put trace=3 to get full details
solgood <- opm(x0, bt.f, bt.g, lower=lo, upper=up, method=meths, control=list(trace=0))
print(solgood)
print(summary(solgood, order=value))
cat("Result for lbfgsb3c\n")
library(lbfgsb3c)
sollbcgood <- lbfgsb3c(x0, bt.f, bt.g, lower=lo, upper=up)
sollbcgood
bady <- bt.badsetup(nn)
lo <- bady$lower
up <- bady$upper
x0 <- bady$x
bady # This has x outside bounds, but method corrects. 
## Should we notify user in a more aggressive fashion?
## optim()
## Put trace=3 to get full details
solbad0 <- opm(x0, bt.f, bt.g, lower=lo, upper=up, method=meths, control=list(trace=0))
print(summary(solbad0, order=value))
sollbcbad <- lbfgsb3c(x0, bt.f, bt.g, lower=lo, upper=up)
sollbcbad
} else {
## Following only works with optimx >= 2019-6.14 since lbfgsb3c not present beforehand
meths <- c("L-BFGS-B", "lbfgsb3c", "lbfgsb3", "Rvmmin", "Rcgmin", "Rtnmin")
nn <- 4
goody <- bt.setup0(nn)
lo <- goody$lower
up <- goody$upper
x0 <- goody$x
goody
## optim()
## Put trace=3 to get full details
solgood <- opm(x0, bt.f, bt.g, lower=lo, upper=up, method=meths, control=list(trace=0))
print(summary(solgood, order=value))
bady <- bt.badsetup(nn)
lo <- bady$lower
up <- bady$upper
x0 <- bady$x
bady # This has x outside bounds, but method corrects. 
## Should we notify user in a more aggressive fashion?
## optim()
## Put trace=3 to get full details
solbad0 <- opm(x0, bt.f, bt.g, lower=lo, upper=up, method=meths, control=list(trace=0))
print(summary(solbad0, order=value))
}


## ---- candlestick-------------------------------------------------------------
# candlestick function
# J C Nash 2011-2-3
cstick.f<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-crossprod(x)
  f<-as.double(r2+alpha/r2)
  return(f)
}

cstick.g<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-as.numeric(crossprod(x))
  g1<-2*x
  g2 <- (-alpha)*2*x/(r2*r2)
  g<-as.double(g1+g2)
  return(g)
}
library(lbfgsb3c)
nn <- 2
x0 <- c(10,10)
lo <- c(1, 1)
up <- c(10,10)
print(x0)
## c2o <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0))
## print(summary(c2o, order=value))
c2l1 <- lbfgsb3c(x0, cstick.f, cstick.g, lower=lo, upper=up)
c2l1

lo <- c(4, 4)
## c2ob <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0))
## print(summary(c2ob, order=value))
c2l2 <- lbfgsb3c(x0, cstick.f, cstick.g, lower=lo, upper=up)
c2l2

## cstick2b <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2b, par.select=1:2, order=value))

## nn <- 100
## x0 <- rep(10, nn)
## up <- rep(10, nn)
## lo <- rep(1e-4, nn)
## cco <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0, kkt=FALSE))
## print(summary(cco, par.select=1:4, order=value))

## ----exrosen------------------------------------------------------------------
# require(funconstrain) ## not in CRAN, so explicit inclusion of this function
# exrosen <- ex_rosen()
# exrosenf <- exrosen$fn
exrosenf <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    fsum <- 0
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        f_p1 <- 10 * (par[p2] - par[p1]^2)
        f_p2 <- 1 - par[p1]
        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
    }
    fsum
}
# exroseng <- exrosen$gr
exroseng <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    grad <- rep(0, n)
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]
        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]
        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
    }
    grad
}

exrosenx0 <- function (n = 20) {
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    rep(c(-1.2, 1), n/2)
}


require(lbfgsb3c)
## require(optimx)

## require(optimx)
for (n in seq(2,12, by=2)) {
  cat("ex_rosen try for n=",n,"\n")
  x0 <- exrosenx0(n)
  lo <- rep(-1.5, n)
  up <- rep(3, n)
  print(x0)
  cat("optim L-BFGS-B\n")
  eo <- optim(x0, exrosenf, exroseng, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=0))
  print(eo)
  cat("lbfgsb3c\n")
  el <- lbfgsb3c(x0, exrosenf, exroseng, lower=lo, upper=up, control=list(trace=0))
  print(el)
##    erfg <- opm(x0, exrosenf, exroseng, method=meths, lower=lo, upper=up)
##    print(summary(erfg, par.select=1:2, order=value))
}

## ---- usingFortran, eval=FALSE------------------------------------------------
#  system("R CMD SHLIB jrosen.f")
#  dyn.load("jrosen.so")
#  is.loaded("rosen")
#  x0 <- as.double(c(-1.2,1))
#  fv <- as.double(-999)
#  n <- as.double(2)
#  testf <- .Fortran("rosen", n=as.integer(n), x=as.double(x0), fval=as.double(fv))
#  testf
#  
#  rrosen <- function(x) {
#    fval <- 0.0
#    for (i in 1:(n-1)) {
#      dx <- x[i + 1] - x[i] * x[i]
#      fval <- fval + 100.0 * dx * dx
#      dx <- 1.0 - x[i]
#      fval <- fval + dx * dx
#    }
#    fval
#  }
#  
#  (rrosen(x0))
#  
#  frosen <- function(x){
#    nn <- length(x)
#    if (nn > 100) { stop("max number of parameters is 100")}
#    fv <- -999.0
#    val <- .Fortran("rosen", n=as.integer(nn), x=as.double(x), fval=as.double(fv))
#    val$fval # NOTE--need ONLY function value returned
#  }
#  # Test the funcion
#  tval <- frosen(x0)
#  str(tval)
#  
#  cat("Run with Nelder-Mead using R function\n")
#  mynm <- optim(x0, rrosen, control=list(trace=0))
#  print(mynm)
#  cat("\n\n Run with Nelder-Mead using Fortran function")
#  mynmf <- optim(x0, frosen, control=list(trace=0))
#  print(mynmf)
#  
#  
#  library(lbfgsb3c)
#  library(microbenchmark)
#  cat("try lbfgsb3c, no Gradient \n")
#  cat("R function\n")
#  tlR<-microbenchmark(myopR <- lbfgsb3c(x0, rrosen, gr=NULL, control=list(trace=0)))
#  print(tlR)
#  print(myopR)
#  cat("Fortran function\n")
#  tlF<-microbenchmark(myop <- lbfgsb3c(x0, frosen, gr=NULL, control=list(trace=0)))
#  print(tlF)
#  print(myop)

