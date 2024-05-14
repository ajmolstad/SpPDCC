## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(tidy = FALSE)

## -----------------------------------------------------------------------------
# install.packages("devtools")
# library(devtools)
# devtools::install_github("ajmolstad/SpPDCC")
library(SpPDCC)

## ----cache=TRUE---------------------------------------------------------------
# ------------------------------   
# Basis covariance matrix  
# ----------------------------- 
p <- 20
Omega <- matrix(0, nrow = p, ncol = p)
for (jj in 1:p) {
  for (kk in 1:p) {
    Omega[jj,kk] <- 0.5^abs(jj-kk) * (abs(jj - kk) < 4)
  }
}
diag(Omega) <- 1
eo <- eigen(Omega)
OmegaSqrt <- eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)

## ----cache=TRUE---------------------------------------------------------------
# -----------------------------  
# Generate log-abundances  
# -----------------------------
set.seed(1)
n <- 100
# create training set
logW <- matrix(rnorm(n*p), nrow = n, ncol = p)%*%OmegaSqrt
W <- exp(logW)
dat <- W/rowSums(W)
dat[1:5, 1:5]

# create validation set
logWval <- matrix(rnorm(n*p), nrow = n, ncol = p)%*%OmegaSqrt
Wval <- exp(logWval)
datval <- Wval/rowSums(Wval)

## ----cache=TRUE, SCCexample---------------------------------------------------
t0 <- SCCpath(dat = dat, datval = datval, 
      nlambda1 = 25, deltalambda1 = 0.01, 
      alpha0 = 10,   tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      silent = FALSE,
      quiet = TRUE,
      path.requires = NULL)
str(t0)
plot(log10(t0$lambda1.vec), t0$valErrs, pch=20, xlab=expression(log[10](lambda[1])), 
     ylab="Validation error")
abline(v = log10(t0$lambda1.vec[which.min(t0$valErrs)]))

## ----cache=TRUE---------------------------------------------------------------
best.ind <- which.min(t0$valErrs)
library(lattice)
levelplot(t0$Omega[1,,,best.ind], col.regions = gray(100:0/100))
levelplot(Omega, col.regions = gray(100:0/100))

## ----cache=TRUE---------------------------------------------------------------
set.seed(1)
t0.cv <- SCCcv(dat, 
      nlambda1 = 25,
      deltalambda1 = 0.01, 
      alpha0 = 10,
      nfolds = 5,   
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = FALSE,
      quiet = TRUE)
str(t0.cv)
plot(log10(t0.cv$lambda1.vec), t0.cv$avg.valErrs, pch=20, xlab=expression(log[10](lambda[1])), ylab="Average CV error")


## ----cache=TRUE---------------------------------------------------------------
# ------------------------------   
# Basis covariance matrix  
# ----------------------------- 
p <- 20
Omega <- array(0, dim=c(2, p, p))
for (jj in 1:p) {
  for (kk in 1:p) {
    Omega[1,jj,kk] <- 0.5^abs(jj-kk) * (abs(jj - kk) < 4)
    Omega[2,jj,kk] <- 0.7^abs(jj-kk) * (abs(jj - kk) < 5)
  }
}

diag(Omega[1,,]) <- 1; diag(Omega[2,,]) <- 1
eo <- eigen(Omega[1,,])
OmegaSqrt1 <- eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
eo <- eigen(Omega[2,,])
OmegaSqrt2 <- eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)

# -----------------------------  
# Generate log-abundances  
# -----------------------------
set.seed(1)
n <- 100
dat <- list(); datval <- list()
# create training set
logW <- matrix(rnorm(n*p), nrow = n, ncol = p)%*%OmegaSqrt1
W <- exp(logW)
dat[[1]] <- W/rowSums(W)
logW <- matrix(rnorm(n*p), nrow = n, ncol = p)%*%OmegaSqrt2
W <- exp(logW)
dat[[2]] <- W/rowSums(W)

# create validation set
logWval <- matrix(rnorm(n*p), nrow = n, ncol = p)%*%OmegaSqrt1
Wval <- exp(logWval)
datval[[1]] <- Wval/rowSums(Wval)
logWval <- matrix(rnorm(n*p), nrow = n, ncol = p)%*%OmegaSqrt2
Wval <- exp(logWval)
datval[[2]] <- Wval/rowSums(Wval)

## ----cache=TRUE---------------------------------------------------------------
t1 <- gSCCpath(dat = dat, 
      datval = datval, 
      weighted = FALSE,
      nlambda1 = 25, nlambda2 = 15, 
      deltalambda1 = 0.001, 
      deltalambda2 = 0.01, 
      alpha0 = 1,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = TRUE,
      quiet = TRUE)
str(t1)
levelplot(t1$valErrs, col.regions = gray(100:0/100))
inds <- which(t1$valErrs == min(t1$valErrs), arr.ind=TRUE)
Omega.est <- t1$Omega[,,,inds[1,1], inds[1,2]]
sum((Omega.est[1,,] - Omega[1,,])^2); sum((Omega.est[2,,] - Omega[2,,])^2)

## ----cache=TRUE---------------------------------------------------------------
t1.cv <- gSCCcv(dat = dat, 
      nfolds = 5, 
      nlambda1 = 20, nlambda2 = 10, 
      deltalambda1 = 0.001, 
      deltalambda2 = 0.01, 
      alpha0 = 1,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = TRUE,
      quiet = TRUE)
str(t1.cv)
sum((t1.cv$Omega.min[1,,] - Omega[1,,])^2); sum((t1.cv$Omega.min[2,,] - Omega[2,,])^2)

## ----cache=TRUE---------------------------------------------------------------
library(lattice)
levelplot(t1.cv$Omega.min[1,,], col.regions = gray(100:0/100))
levelplot(Omega[1,,], col.regions = gray(100:0/100))

## ----cache=TRUE---------------------------------------------------------------
t1w.path <- gSCCpath(dat = dat, datval = datval,
      nlambda1 = 20, nlambda2 = 10, 
      deltalambda1 = 0.001, 
      deltalambda2 = 0.01, 
      weighted = TRUE,
      alpha0 = 1,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = TRUE,
      quiet = TRUE)
str(t1w.path)
sum((t1.cv$Omega.min[1,,] - Omega[1,,])^2); sum((t1.cv$Omega.min[2,,] - Omega[2,,])^2)

