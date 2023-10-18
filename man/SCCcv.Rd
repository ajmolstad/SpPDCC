\name{SCCcv}
\alias{SCCcv}
\title{SCCcv}
\description{Fit the solution path and perform cross-validation for the sparse and positive definite estimator of multiple basis covariance matrices from compositional data.}
\usage{
SCCcv(dat, nlambda1 = 25,
      deltalambda1 = 0.01, nfolds = 5,  alpha0 = 10,  
      tau = 1e-4, max.iter = 5e3, tol.obj = 1e-6, tol.diff = 1e-6,
      tol.obj.con = 1e-6, tol.diff.con = 1e-6, silent = FALSE, quiet = TRUE)
}



\arguments{
\item{dat}{Training data: a matrix of dimension \eqn{n \times p}. Each row of the matrix must have positive entries which sum to one.}
\item{nlambda1}{Number of candidate tuning parameters for the lasso-type penalty applied to each basis covariance matrix. }
\item{deltalambda1}{The ratio of minimum to maximum \eqn{\lambda_1} tuning parameter values: must be between zero and one. Decrease if validation errors are minimized for smallest candidate \eqn{\lambda_1} value.}
\item{nfolds}{Positive integer indicating the number of folds to be used for cross-validation.}
For all others, see the example document on \url{https://github.com/ajmolstad/SpPDCC}. 
}

\value{
\item{\code{Omega}}{An array of dimension \eqn{1 \times p \times p \times \texttt{nlambda1}}. The first three dimensions correspond to the estimated covariance matrix for the different candidate tuning parameter values.}
\item{\code{valErrs}}{A matrix of dimension \eqn{\texttt{nfolds} \times \texttt{nlambda1}} with validation errors for each candidate tuning parameter.}
\item{\code{lambda1.vec}}{The candidate tuning parameters for \eqn{\lambda_1}. }
\item{\code{avg.valErrs}}{A vector of length \eqn{\texttt{nlambda1}} with the cross-validation errors averaged over the folds.}
\item{\code{Omega.min}}{An array of dimension \eqn{1 \times p \times p} with the covariance matrix estimate using the tuning parameter which minimized the averaged cross-validation error. }
\item{\code{fold.id}}{The fold IDs used for cross-validation.}
}

