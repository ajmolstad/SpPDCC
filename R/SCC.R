# -----------------------------------------------------
# Auxillary functions
# ----------------------------------------------------- 
remove.diags <- function(X) {
  out <- X
  for (j in 1:dim(X)[1]) {
    diag(out[j,,]) <- 0
  }
  return(out)
}

obj.func <- function(Ts, Omega, lambda1, lambda2) {
  out <- 0
  for (j in 1:dim(Ts)[1]) {
    out <- out + sum((Ts[j,,] - matrix(diag(Omega[j,,]), nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=TRUE)  - matrix(diag(Omega[j,,]),nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=FALSE) + 2*Omega[j,,])^2)
  }
  return(out + lambda1*sum(abs(remove.diags(Omega))) + lambda2*sum(apply(remove.diags(Omega), c(2,3), function(x){sqrt(sum(x^2))})))
}


loss.func <- function(Ts, Omega) {
  out <- 0
  for (j in 1:dim(Ts)[1]) {
    out <- out + sum((Ts[j,,] - matrix(diag(Omega[j,,]), nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=TRUE)  - matrix(diag(Omega[j,,]),nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=FALSE) + 2*Omega[j,,])^2)
  }
  return(out)
}


prox.function <- function(Y, lambda1, lambda2) {
  if(lambda2!=0){
    temp <- pmax(abs(Y) - lambda1, 0)*sign(Y)
    out <- array(0, dim = dim(Y))
    for (kk in 1:(dim(Y)[2]-1)) {
      for (jj in (kk+1):dim(Y)[3]) {
        if (kk != jj) {
          t0 <- sqrt(sum(temp[,kk,jj]^2))
          out[,jj,kk] <- max(1 - lambda2/t0, 0)*temp[,kk,jj]
        } 
      }
    }
    for (j in 1:dim(Y)[1]) {
      out[j,,] <- out[j,,] + t(out[j,,])
      diag(out[j,,]) <- diag(Y[j,,])
    }
  } else {
    out <- pmax(abs(Y) - lambda1, 0)*sign(Y)
    for (j in 1:dim(Y)[1]) {
      diag(out[j,,]) <- diag(Y[j,,])
    } 
  }
  return(out)
}



grad.loss <- function(Ts, Omega) {
  out <- array(0, dim=dim(Ts))
  for (j in 1:dim(Ts)[1]) {
    for (k in 1:dim(Ts)[2]) {
      for (l in k:dim(Ts)[3]) {
        if (k == l) {
          for (mm in c(1:dim(Ts)[2])[-k]) {
            out[j,k,l] <- out[j,k,l] + 4*Omega[j,k,k] - 4*(Ts[j,k,mm] - Omega[j,mm,mm] + 2*Omega[j,k,mm])
          }
        } else {
          out[j,k,l] <- 4*(Ts[j,k,l] - Omega[j,k,k] - Omega[j,l,l]) + 8*Omega[j,k,l]
        }  
      }
    }
    out[j,,] <- out[j,,] + t(out[j,,])
    diag(out[j,,]) <- diag(out[j,,])/2
  }
  return(out)
}


# --------------------------------------------------------------
# Prox-prox gradient descent algorithm
# Used when unconstrained estimator does not satsify PDness
# --------------------------------------------------------------

SCC.constrained <- function(Ts,
  X.init = NULL,
  Z.init = NULL,
  U.init = NULL,
  alpha0 = 1, 
  tau = 1e-4,
  max.iter = 5e3,
  lambda1,
  lambda2,
  tol.obj = 1e-8,
  tol.diff = 1e-8,
  quiet = TRUE) {

  alpha <- alpha0
  X.new <- X.init
  Z.new <- Z.init
  U.new <- U.init
  loss.prev <- loss.func(Ts, Z.new)
  obj.func <- rep(0, max.iter)

  for (iteration in 1:max.iter) {
    
    grad.temp <- grad.loss(Ts, Z.new)
    lineSearch <- TRUE

    while (lineSearch) {
      X.temp <- Z.new - alpha*U.new - alpha*grad.temp
      X.new <- prox.function(X.temp, alpha*lambda1, alpha*lambda2)
      loss.temp <- loss.func(Ts, X.new)
      if (loss.temp <= loss.prev + sum(grad.temp*(X.new - Z.new)) + sum((X.new - Z.new)^2)/(2*alpha)) {
        lineSearch <- FALSE
      } else {
        alpha <- alpha*.5
      }
    }
    
    
    Z.new <- X.new + alpha*U.new
    for (kk in 1:dim(Z.new)[1]) {
      eo <- eigen(Z.new[kk,,])
      if (!all(eo$val > tau)) {
        Z.new[kk,,] <- crossprod(t(eo$vec), pmax(eo$val, tau)*t(eo$vec))
        #no.Need <- FALSE
      } 
    }

    U.temp <- array(0, dim=dim(Z.new))
    for (kk in 1:dim(Z.new)[1]) {
      U.temp[kk,,] <- U.new[kk,,] + (X.new[kk,,] - Z.new[kk,,])/alpha
      U.temp[kk,,] <- (U.temp[kk,,] + t(U.temp[kk,,]))/2
    }

    U.new <- U.temp
    loss.prev <- loss.func(Ts, Z.new)

    obj.func[iteration] <- loss.prev + lambda1*sum(abs(remove.diags(Z.new))) + lambda2*sum(apply(remove.diags(Z.new), c(2,3), function(x){sqrt(sum(x^2))}))
    if (!quiet) {
      cat(obj.func[iteration], "\n")
    }
    if (iteration > 5) {
      if (all(abs(obj.func[iteration:(iteration-2)] - obj.func[(iteration-1):(iteration-3)]) < tol.obj*abs(obj.func[iteration]))) {
        if (max(abs(X.new - Z.new)) < tol.diff) {
          break
        }
      }
    }
  }

  return(list("U.new" = U.new, "Z.new" = Z.new, "X.new" = X.new))
}


# -------------------------------------------------------------------
# Acc prox gradient descent algorithm for unconstrained estimator
# ------------------------------------------------------------------
SCC <- function(Ts,
  X.init = NULL,
  alpha0 = 1, 
  tau = 1e-4,
  max.iter = 5e3,
  lambda1,
  lambda2,
  tol.obj = 1e-8,
  tol.diff = 1e-8,
  quiet = TRUE) {

  alpha <- alpha0
  X.current <- X.init
  X.prev <- X.init
  loss.prev <- loss.func(Ts, X.prev)
  obj.func <- rep(0, max.iter)

  for (iteration in 1:max.iter) {
    
    X.search <- X.current + (iteration/(iteration + 3))*(X.current - X.prev)
    grad.temp <- grad.loss(Ts, X.search)
    lineSearch <- TRUE

    while (lineSearch) {

      X.temp <- X.search - alpha*grad.temp
      X.new <- prox.function(X.temp, alpha*lambda1, alpha*lambda2)
      loss.temp <- loss.func(Ts, X.new)
      loss.search <- loss.func(Ts, X.search)
      if (loss.temp <= loss.search + sum(grad.temp*(X.new - X.search)) + sum((X.new - X.search)^2)/(2*alpha)) {
        lineSearch <- FALSE
        X.prev <- X.current
        X.current <- X.new
        loss.prev <- loss.temp
      } else {
        alpha <- alpha*.5
      }
    }
    

    obj.func[iteration] <- loss.prev + lambda1*sum(abs(remove.diags(X.current))) + lambda2*sum(apply(remove.diags(X.current), c(2,3), function(x){sqrt(sum(x^2))}))
    if (!quiet) {
      cat(obj.func[iteration], "\n")
    }
    if (iteration > 5) {
      if (all(abs(obj.func[iteration:(iteration-2)] - obj.func[(iteration-1):(iteration-3)]) < tol.obj*abs(obj.func[iteration]))) {
        if (max(abs(X.current - X.prev)) < tol.diff) {
          break
        }
      }
    }
  }

  return(list("X.current" = X.current))
}


# -----------------------------------------------------
# Compute solution path and compare to validation 
# -----------------------------------------------------
# Args: 
#     - Ts: sample compositional covariance (list, 1 x p x p)
#     - TsVal: validation set compositional covariance (list, 1 x p x p)
#     - nlambda1: number of candidate lambda1
#     - deltalambda1: ratio of lambda1_max to lambda1_min
#     - alpha0: starting step size (default 10)
#     - tau: lower bound on eigenvalues of solution
#     - max.iter: maximum AccPGD or ProxProxGD iterations
#     - tol.obj: % change in objective function for convergence
#     - tol.diff: max-abs change in iterates for convergence
#     - silent: print validation loss? 
#     - quiet: print objective function after each iterate?
#     - path.requires: how many tuning parameters must be tried 
#         we terminate path when validation loss increases for five
#         successive candidate tuning parameters
# -----------------------------------------------------

SCCpath <- function(dat, 
      datval = NULL, 
      nlambda1 = 25,
      deltalambda1 = 0.01, 
      lambda1.vec = NULL, 
      alpha0 = 10,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = TRUE,
      quiet = TRUE,
      path.requires = NULL){

  # -------------------------------------------
  # Preliminary checks 
  # -------------------------------------------
  if (is.null(path.requires)) {
    path.requires <- nlambda1
  }
  if (!is.matrix(dat)) {
    stop("dat needs to an n-by-p matrix be a matrix with p-dimensional compositional observations organized by row")
  }
  if (is.null(datval)) {
    datval <- dat
    returnVals <- FALSE
  } else {
    returnVals <- TRUE
  }

  if (silent) {
    quiet <- TRUE
  }
  

  p <- dim(dat)[2]
  J <- 1

  # -------------------------------------------
  # Compute sample variation matrices
  # -------------------------------------------
  Ts <- array(0, dim=c(1, p, p))
  Tsval <- array(0, dim=c(1, p, p))
  n <- dim(dat)[1]
  nval <- dim(datval)[1]

  for (kk in 1:J) {
    # ------------------------------
    # training set 
    # ------------------------------
    gMat <- matrix(0, p, p)
    for (jj in 1:p) {
      for (ll in 1:p) {
         gMat[jj,ll] <- sum(log(dat[,jj]/dat[,ll]))/n
      }
    }
    for (jj in 1:p) {
      for (ll in 1:p) {
        Ts[kk,jj,ll] <- sum((log(dat[,jj]/dat[,ll]) - gMat[jj,ll])^2)/n
      }
    }
    # ------------------------------
    # validation set 
    # ------------------------------
    gMat <- matrix(0, p, p)
    for (jj in 1:p) {
      for (ll in 1:p) {
         gMat[jj,ll] <- sum(log(datval[,jj]/datval[,ll]))/nval
      }
    }
    for (jj in 1:p) {
      for (ll in 1:p) {
        Tsval[kk,jj,ll] <- sum((log(datval[,jj]/datval[,ll]) - gMat[jj,ll])^2)/nval
      }
    }
  }

  # -------------------------------------------
  # Algorithm to compute diagonal estimator
  # -------------------------------------------
  getDiagT <- function(T, max.iter = 1e4){
    p <- dim(T)[1]
    omega <- rep(1, dim(T)[1])
    omega.prev <- omega
    for (kk in 1:max.iter) {
      for (jj in 1:p) {
        omega[jj] <- sum(T[jj,-jj] - omega[-jj])/(p-1)
      }
      if (sum(abs(omega - omega.prev)) < 1e-12) {
        break
      } else {
        omega.prev <- omega
      }
    }
    return(omega)
  }

  getDiagTs <- function(Ts){
    Omega <- array(0, dim(Ts))
    for (kk in 1:dim(Ts)[1]) {
      diag(Omega[kk,,]) <- getDiagT(Ts[kk,,])
    }
    return(Omega)
  }

  # ---------------------------------------------
  # Find candidate tuning parameters
  # ---------------------------------------------
  out <- getDiagTs(Ts)
  t0 <- grad.loss(Ts, out)
  if(is.null(lambda1.vec)){
    lambda1max.candidate <- max(remove.diags(abs(t0)))
    lambda1.vec <- 10^seq(log10(lambda1max.candidate), log10(deltalambda1*lambda1max.candidate), length=nlambda1)
  } 
  nlambda1 <- length(lambda1.vec)
  U.init <- array(0, dim= dim(Ts))
  X.init <- Ts
  Z.init <- Ts
  Omega.mat <- array(0, dim=c(dim(Ts)[1], dim(Ts)[2], dim(Ts)[3], nlambda1))
  results <- rep(Inf, nlambda1)
  unconstrained <- TRUE

  # ----------------------------------------------
  # Compute solution path
  # ----------------------------------------------
  for (jj in 1:nlambda1) {
    if (unconstrained) {
      out <- SCC(Ts,
        X.init = X.init,
        alpha0 = alpha0, 
        tau = tau,
        max.iter = max.iter,
        lambda1 = lambda1.vec[jj],
        lambda2 = 0,
        tol.obj = tol.obj,
        tol.diff = tol.diff,
        quiet = quiet)
      if (min(eigen(out$X.current[1,,])$val) >= tau) {
        Omega.mat[,,,jj] <- out$X.current[1,,]
        X.init <- out$X.current
        Z.init <- out$X.current
        results[jj] <- loss.func(Tsval, out$X.current)
      } else {
        unconstrained <- FALSE
      }
    }

    if (!unconstrained) {
      out <- SCC.constrained(Ts,
          X.init = X.init,
          Z.init = Z.init,
          U.init = U.init,
          alpha0 = alpha0, 
          tau = tau,
          max.iter = max.iter,
          lambda1 = lambda1.vec[jj],
          lambda2 = 0,
          tol.obj = tol.obj.con,
          tol.diff = tol.diff.con,
          quiet = quiet)

      X.init <- out$X.new
      Z.init <- out$Z.new
      U.init <- out$U.new
      results[jj] <- loss.func(Tsval, out$Z.new)
      Omega.mat[,,,jj] <- out$X.new
    }

    if (jj > path.requires) {
      if (all(results[jj:(jj-4)] > results[(jj-1):(jj-5)])) {
        break
      }
    }
    unconstrained <- TRUE 
    if (!silent) {
      cat(results[jj], "\n")
    }
  }

  if(returnVals){
    return(list("Omega" = Omega.mat, 
                "lambda1.vec" = lambda1.vec,
                "valErrs" = results))
  } else {
    return(list("Omega" = Omega.mat, 
                "lambda1.vec" = lambda1.vec))
  }

}



SCCcv <- function(dat, 
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
      quiet = TRUE){


  # -------------------------------------
  # Preliminary checks
  # -------------------------------------
  if (!is.matrix(dat)) {
    stop("dat needs to an n-by-p matrix be a matrix with p-dimensional compositional observations organized by row")
  }
  if (silent) {
    quiet <- TRUE
  }

  n <- dim(dat)[1]
  p <- dim(dat)[2]
  J <- 1

  # -----------------------------------
  # Fit model to full dataset 
  # ------------------------------------
  out <- SCCpath(dat = dat, 
      datval = NULL,
      nlambda1 = nlambda1,
      deltalambda1 = deltalambda1, 
      alpha0 = alpha0,  
      tau = tau,
      max.iter = max.iter,
      tol.obj = tol.obj,
      tol.diff = tol.diff,
      tol.obj.con = tol.obj.con,
      tol.diff.con = tol.diff.con,
      silent = TRUE,
      quiet = TRUE)

  lambda1.vec <- out$lambda1
  if (!silent) {
    cat("Beginning cross-validation; printing left-out fold errors....", "\n")
  }

  # -------------------------------------------
  # Partition the data into nfolds 
  # -------------------------------------------
  fold.id <- sample(rep(1:nfolds, length = n))
  valErrs <- matrix(0, nrow = nfolds, ncol = nlambda1)
  fold.proportion <- table(fold.id)/n
  # -------------------------------------------
  # Perform cross-validation 
  # -------------------------------------------
  for(cv.fold in 1:nfolds){
    dat.train <- dat[-which(fold.id == cv.fold), ]
    dat.val <- dat[which(fold.id == cv.fold), ]
    cv.out <- SCCpath(dat = dat.train, datval = dat.val, 
      lambda1.vec = lambda1.vec,
      alpha0 = alpha0,  
      tau = tau,
      max.iter = max.iter,
      tol.obj = tol.obj,
      tol.diff = tol.diff,
      tol.obj.con = tol.obj.con,
      tol.diff.con = tol.diff.con,
      silent = silent,
      quiet = quiet)
    valErrs[cv.fold,] <- cv.out$valErrs
    if (!silent) {
      cat("# ----------------------------------", "\n")
      cat("Through fold ", cv.fold, " of ", nfolds, "\n")
      cat("# ----------------------------------", "\n")
    }
  }

  valErrs.avg <- colSums(matrix(fold.proportion, nrow=nfolds, ncol=nlambda1, byrow=FALSE)*valErrs)

  return(list("Omega" = out$Omega, 
              "valErrs" = valErrs, 
              "avg.valErrs" = valErrs.avg,
              "Omega.min" = out$Omega[1,,,which.min(valErrs.avg)],
              "fold.id" = fold.id,
              "lambda1.vec" = lambda1.vec
              ))
}


