
remove.diags <- function(X) {
  out <- X
  for (j in 1:dim(X)[1]) {
    diag(out[j,,]) <- 0
  }
  return(out)
}

obj.func.w <- function(Ts, ns, Omega, lambda1, lambda2){
  out <- 0
  N <- sum(ns)
  for(j in 1:dim(Ts)[1]){
    out <- out + (ns[j]/N)*sum((Ts[j,,] - matrix(diag(Omega[j,,]), nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=TRUE)  - matrix(diag(Omega[j,,]),nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=FALSE) + 2*Omega[j,,])^2)
  }
  return(out + lambda1*sum(abs(remove.diags(Omega))) + lambda2*sum(apply(remove.diags(Omega), c(2,3), function(x){sqrt(sum(x^2))})))
}


loss.func.w <- function(Ts, ns, Omega){
  out <- 0
  N <- sum(ns)
  for(j in 1:dim(Ts)[1]){
    out <- out + (ns[j]/N)*sum((Ts[j,,] - matrix(diag(Omega[j,,]), nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=TRUE)  - matrix(diag(Omega[j,,]),nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=FALSE) + 2*Omega[j,,])^2)
  }
  return(out)
}


grad.loss.w <- function(Ts, ns, Omega){
  out <- array(0, dim=dim(Ts))
  N <- sum(ns)
  for(j in 1:dim(Ts)[1]){
    for(k in 1:dim(Ts)[2]){
      for(l in k:dim(Ts)[3]){
        if(k == l){
          for(mm in c(1:dim(Ts)[2])[-k]){
            out[j,k,l] <- out[j,k,l] + 4*Omega[j,k,k] - 4*(Ts[j,k,mm] - Omega[j,mm,mm] + 2*Omega[j,k,mm])
          }
        } else {
          out[j,k,l] <- 4*(Ts[j,k,l] - Omega[j,k,k] - Omega[j,l,l]) + 8*Omega[j,k,l]
        }  
      }
    }
    out[j,,] <- out[j,,] + t(out[j,,])
    diag(out[j,,]) <- diag(out[j,,])/2
    out[j,,] <- (ns[j]/N)*out[j,,]
  }
  return(out)
}



prox.function.array <- function(Y, lambda1, lambda2) {
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
  return(out)
}




gSCCw.constrained <- function(
  Ts, 
  ns,
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
  quiet = TRUE
) {

  alpha <- alpha0
  X.new <- X.init
  Z.new <- Z.init
  U.new <- U.init
  loss.prev <- loss.func.w(Ts, ns, Z.new)
  obj.func <- rep(0, max.iter)

  for(iteration in 1:max.iter){
    
    grad.temp <- grad.loss.w(Ts, ns, Z.new)
    lineSearch <- TRUE

    while(lineSearch){
      X.temp <- Z.new - alpha*U.new - alpha*grad.temp
      X.new <- prox.function(X.temp, alpha*lambda1, alpha*lambda2)
      loss.temp <- loss.func.w(Ts, ns, X.new)
      if(loss.temp <= loss.prev + sum(grad.temp*(X.new - Z.new)) + sum((X.new - Z.new)^2)/(2*alpha)){
        lineSearch <- FALSE
      } else {
        alpha <- alpha*.5
      }
    }
    
    
    Z.new <- X.new + alpha*U.new
    for(kk in 1:dim(Z.new)[1]){
      eo <- eigen(Z.new[kk,,])
      if(!all(eo$val > tau)){
        Z.new[kk,,] <- crossprod(t(eo$vec), pmax(eo$val, tau)*t(eo$vec))
      } 
    }

    U.temp <- array(0, dim=dim(Z.new))
    for(kk in 1:dim(Z.new)[1]){
      U.temp[kk,,] <- U.new[kk,,] + (X.new[kk,,] - Z.new[kk,,])/alpha
      U.temp[kk,,] <- (U.temp[kk,,] + t(U.temp[kk,,]))/2
    }

    U.new <- U.temp
    loss.prev <- loss.func.w(Ts, ns, Z.new)

    obj.func[iteration] <- loss.prev + lambda1*sum(abs(remove.diags(Z.new))) + lambda2*sum(apply(remove.diags(Z.new), c(2,3), function(x){sqrt(sum(x^2))}))
    if(!quiet){
      cat(obj.func[iteration], "\n")
    }
    if(iteration > 5){
      if(all(abs(obj.func[iteration:(iteration-2)] - obj.func[(iteration-1):(iteration-3)]) < tol.obj*abs(obj.func[iteration]))){
        if(max(abs(X.new - Z.new)) < tol.diff){
          break
        }
      }
    }
  }

  return(list("U.new" = U.new, "Z.new" = Z.new, "X.new" = X.new))
}



gSCCw <- function(
  Ts, 
  ns,
  X.init = NULL,
  alpha0 = 1, 
  tau = 1e-4,
  max.iter = 5e3,
  lambda1,
  lambda2,
  tol.obj = 1e-8,
  tol.diff = 1e-8,
  quiet = TRUE
) {

  alpha <- alpha0
  X.current <- X.init
  X.prev <- X.init
  loss.prev <- loss.func.w(Ts, ns, X.prev)
  obj.func <- rep(0, max.iter)

  for(iteration in 1:max.iter){
    
    X.search <- X.current + (iteration/(iteration + 3))*(X.current - X.prev)
    grad.temp <- grad.loss.w(Ts, ns, X.search)
    lineSearch <- TRUE

    while (lineSearch) {

      X.temp <- X.search - alpha*grad.temp
      X.new <- prox.function(X.temp, alpha*lambda1, alpha*lambda2)
      loss.temp <- loss.func.w(Ts, ns, X.new)
      loss.search <- loss.func.w(Ts, ns, X.search)
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



gSCCwpath <- function(dat, 
      datval = NULL, 
      nlambda1 = 25, nlambda2 = 25, 
      deltalambda1 = 0.01, 
      deltalambda2 = 0.01, 
      lambda1.mat = NULL, 
      lambda2.vec = NULL,
      alpha0 = 10,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = FALSE,
      quiet = TRUE,
      path.requires = NULL) {

  # -------------------------------------------
  # Preliminary checks 
  # -------------------------------------------
  if (!is.list(dat)) {
    stop("dat needs to a list containing n_h-by-p matrices with p-dimensional compositional observations organized by row")
  }
  if (is.null(datval)) {
    datval <- dat
    return.valErrs <- FALSE
  } else {
    return.valErrs <- TRUE
  }
  if (silent) {
    quiet <- TRUE
  }
  if (is.null(path.requires)) {
    path.requires <- Inf
  }

  p <- dim(dat[[1]])[2]
  J <- length(dat)

  # -------------------------------------------
  # Compute sample variation matrices
  # -------------------------------------------
  Ts <- array(0, dim=c(J, p, p))
  Tsval <- array(0, dim=c(J, p, p))
  ns <- rep(0, J)
  nsval <- rep(0, J)
  for (kk in 1:J) {

    ns[kk] <- dim(dat[[kk]])[1]
    nsval[kk] <- dim(datval[[kk]])[1]

    gMat <- matrix(0, p, p)
    for (jj in 1:p) {
        for (ll in 1:p) {
           gMat[jj,ll] <- sum(log(dat[[kk]][,jj]/dat[[kk]][,ll]))/dim(dat[[kk]])[1]
        }
    }
    for (jj in 1:p) {
        for (ll in 1:p) {
          Ts[kk,jj,ll] <- sum((log(dat[[kk]][,jj]/dat[[kk]][,ll]) - gMat[jj,ll])^2)/dim(dat[[kk]])[1]
        }
    }

    gMat <- matrix(0, p, p)
    for (jj in 1:p) {
        for (ll in 1:p) {
           gMat[jj,ll] <- sum(log(datval[[kk]][,jj]/datval[[kk]][,ll]))/dim(datval[[kk]])[1]
        }
    }
    for (jj in 1:p) {
        for (ll in 1:p) {
          Tsval[kk,jj,ll] <- sum((log(datval[[kk]][,jj]/datval[[kk]][,ll]) - gMat[jj,ll])^2)/dim(datval[[kk]])[1]
        }
    }

  }

  # -------------------------------------------------------------------
  # Algorithm to compute estimator when all off-diagonals are zero
  # --------------------------------------------------------------------
  getDiagT <- function(T, max.iter = 1e4) {
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

  getDiagTs <- function(Ts) {
    Omega <- array(0, dim(Ts))
    for (kk in 1:dim(Ts)[1]) {
      diag(Omega[kk,,]) <- getDiagT(Ts[kk,,])
    }
    return(Omega)
  }


  # ------------------------------------------------------------------
  # Get candidate tuning parameter set
  # ------------------------------------------------------------------
  if (is.null(lambda1.mat) | is.null(lambda2.vec)) {

    out <- getDiagTs(Ts)
    t0 <- grad.loss.w(Ts, ns, out)
    lambda2max <- max(apply(t0, c(2,3), function(x){sqrt(sum(x^2))}))
    lambda2.vec <- 10^seq(log10(lambda2max), log10(deltalambda2*lambda2max), length=nlambda2+2)[-c(1,2)]
    lambda1.mat <- matrix(0, nrow=nlambda2, ncol=nlambda1)
    lambda1max.candidate <- max(remove.diags(abs(t0)))
    lambda1.candidates <- 10^seq(log10(lambda1max.candidate), log10(.01*lambda1max.candidate), length=100)
    U.init <- array(0, dim= dim(Ts))
    X.init <- Ts
    Z.init <- Ts
    unconstrained <- TRUE

    for (kk in 1:nlambda2) {

      for (jj in 1:100) {

        if (unconstrained) {

         out <- gSCCw(Ts, ns,
          X.init = X.init,
          alpha0 = alpha0, 
          tau = tau,
          max.iter = max.iter,
          lambda1 = lambda1.candidates[jj],
          lambda2 = lambda2.vec[kk],
          tol.obj = tol.obj,
          tol.diff = tol.diff,
          quiet = quiet)

          for(ll in 1:dim(out$X.current)[1]){
            if (min(eigen(out$X.current[ll,,])$val) < tau) {
                unconstrained <- FALSE
            }
          }
          X.init <- out$X.current
          Z.init <- out$X.current
        }

        if (!unconstrained) {
          out <- gSCCw.constrained(Ts, ns, 
              X.init = X.init,
              Z.init = Z.init,
              U.init = U.init,
              alpha0 = alpha0, 
              tau = tau,
              max.iter = max.iter,
              lambda1 = lambda1.candidates[jj],
              lambda2 = lambda2.vec[kk],
              tol.obj = tol.obj.con,
              tol.diff = tol.diff.con,
              quiet = quiet)

          X.init <- out$X.new
          Z.init <- out$Z.new
          U.init <- out$U.new
          Omega.mat[,,,jj] <- out$X.new
        }

        if (sum(remove.diags(X.init)!=0) > 0) {
          break
        } 
      }
      if (jj > 1) {
        lambda1.mat[kk,] <- 10^seq(log10(lambda1.candidates[jj-1]), log10(deltalambda1*lambda1.candidates[jj-1]), length=nlambda1)
      } else {
        lambda1.mat[kk,] <- 10^seq(log10(lambda1.candidates[jj]), log10(deltalambda1*lambda1.candidates[jj]), length=nlambda1)
      }
    }
  } else {
    nlambda1 <- dim(lambda1.mat)[2]
    nlambda2 <- dim(lambda1.mat)[1]
  }

  # --------------------------------------------------
  # Compute solution path
  # --------------------------------------------------
  Omega.mat <- array(0, dim=c(dim(Ts)[1], dim(Ts)[2], dim(Ts)[3], nlambda1, nlambda2))
  results <- matrix(Inf, nrow=nlambda1, ncol=nlambda2)
  unconstrained <- TRUE
  U.init <- array(0, dim= dim(Ts))
  X.init <- Ts
  Z.init <- Ts
  
  for (kk in 1:nlambda2) {
    for (jj in 1:nlambda1) {

      if (unconstrained) {
       out <- gSCCw(Ts, ns,
        X.init = X.init,
        alpha0 = alpha0, 
        tau = tau,
        max.iter = max.iter,
        lambda1 = lambda1.mat[kk,jj],
        lambda2 = lambda2.vec[kk],
        tol.obj = tol.obj,
        tol.diff = tol.diff,
        quiet = quiet)

        for (ll in 1:dim(out$X.current)[1]) {
          if (min(eigen(out$X.current[ll,,])$val) < tau) {
              unconstrained <- FALSE
          }
        }
        X.init <- out$X.current
        Z.init <- out$X.current
      }

      if (!unconstrained) {
        out <- gSCCw.constrained(Ts, ns,
            X.init = X.init,
            Z.init = Z.init,
            U.init = U.init,
            alpha0 = alpha0, 
            tau = tau,
            max.iter = max.iter,
            lambda1 = lambda1.mat[kk,jj],
            lambda2 = lambda2.vec[kk],
            tol.obj = tol.obj.con,
            tol.diff = tol.diff.con,
            quiet = quiet)
        X.init <- out$X.new
        Z.init <- out$Z.new
        U.init <- out$U.new
      }

      Omega.mat[,,,jj,kk] <- X.init
      results[jj,kk] <- loss.func.w(Tsval, nsval, Omega.mat[,,,jj,kk])

      if (jj > path.requires) {
        if (all(results[jj:(jj-3),kk] > results[(jj-1):(jj-4),kk])) {
          break
        }
      }

      if (!silent) {
        cat (results[jj,kk], "\n")
      }

    }
  }

  if(return.valErrs){
    return(list("Omega" = Omega.mat, 
                "valErrs" = results,
                "lambda1.mat" = lambda1.mat,
                "lambda2.vec" = lambda2.vec))
  } else {
    return(list("Omega" = Omega.mat, 
                "lambda1.mat" = lambda1.mat,
                "lambda2.vec" = lambda2.vec))
  }
}



gSCCwcv <- function(dat, 
      nlambda1 = 25, 
      nlambda2 = 25, 
      deltalambda1 = 0.01, 
      deltalambda2 = 0.01, 
      nfolds = 5,
      alpha0 = 10,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-4,
      tol.diff = 1e-4,
      tol.obj.con = 1e-6,
      tol.diff.con = 1e-6,
      silent = FALSE,
      quiet = TRUE){

  # -------------------------------------
  # Preliminary checks
  # -------------------------------------
  if (!is.list(dat)) {
    stop("dat needs to a list consisting of n_h-by-p matrices with p-dimensional compositional observations organized by row")
  }
  if (silent) {
    quiet <- TRUE
  }
  J <- length(dat)

  # -----------------------------------
  # Fit model to full dataset 
  # ------------------------------------
  out <- gSCCwpath(dat = dat, 
      datval = NULL,
      nlambda1 = nlambda1,
      nlambda2 = nlambda2, 
      deltalambda1 = deltalambda1, 
      deltalambda2 = deltalambda2, 
      alpha0 = alpha0,  
      tau = tau,
      max.iter = max.iter,
      tol.obj = tol.obj,
      tol.diff = tol.diff,
      tol.obj.con = tol.obj.con,
      tol.diff.con = tol.diff.con,
      silent = TRUE,
      quiet = TRUE)

  lambda1.mat <- out$lambda1
  lambda2.vec <- out$lambda2
  if (!silent) {
    cat("Beginning cross-validation; printing left-out fold errors....", "\n")
  }

  # -------------------------------------------
  # Partition the data into nfolds 
  # -------------------------------------------
  fold.ids <- list()
  fold.proportion <- rep(0, nfolds)
  for (kk in 1:J){
    fold.ids[[kk]] <- sample(rep(1:nfolds, length = dim(dat[[kk]])[1]))
  }
  # for(kk in 1:nfolds){
  #   fold.proportion[kk] <- sum(unlist(lapply(fold.ids, function(x){length(which(x == kk))})))
  # }
  # fold.proportion <- fold.proportion/sum(fold.proportion)
  valErrs <- array(0, dim = c(nfolds, nlambda1, nlambda2))

  # -------------------------------------------
  # Perform cross-validation 
  # -------------------------------------------
  for(cv.fold in 1:nfolds){
    dat.train <- list(); dat.val <- list()
    for(kk in 1:J){
      dat.train[[kk]] <- dat[[kk]][-which(fold.ids[[kk]] == cv.fold), ]
      dat.val[[kk]] <- dat[[kk]][which(fold.ids[[kk]] == cv.fold), ]
    }
    cv.out <- gSCCwpath(dat = dat.train, 
      datval = dat.val, 
      lambda1.mat = lambda1.mat,
      lambda2.vec = lambda2.vec,
      alpha0 = alpha0,  
      tau = tau,
      max.iter = max.iter,
      tol.obj = tol.obj,
      tol.diff = tol.diff,
      tol.obj.con = tol.obj.con,
      tol.diff.con = tol.diff.con,
      silent = silent,
      quiet = quiet)
    valErrs[cv.fold,,] <- cv.out$valErrs
    if (!silent) {
      cat("# ----------------------------------", "\n")
      cat("Through fold ", cv.fold, " of ", nfolds, "\n")
      cat("# ----------------------------------", "\n")
    }
  }

  valErrs.avg <- apply(valErrs, c(2,3), sum)
  inds <- which(valErrs.avg == min(valErrs.avg), arr.ind=TRUE)

  return(list("Omega" = out$Omega, 
              "valErrs" = valErrs, 
              "avg.valErrs" = valErrs.avg,
              "fold.ids" = fold.ids,
              "lambda1.mat" = lambda1.mat,
              "lambda2.vec" = lambda2.vec,
              "Omega.min" = out$Omega[,,,inds[1,1], inds[1,2]]))
}

