\name{SCCcv}
\alias{SCCcv}
\title{Cross-validation function for fitting and tuning estimator sparse and positive definite basis covariance matrix}
\description{Fit the solution path and perform cross-validation for the sparse and positive definite estimator of multiple basis covariance matrices from compositional data.}
\usage{
SCCcv(dat, nlambda1 = 25,
      deltalambda1 = 0.01, alpha0 = 10,  nfolds = 5, 
      tau = 1e-4, max.iter = 5e3, tol.obj = 1e-6, tol.diff = 1e-6,
      tol.obj.con = 1e-6, tol.diff.con = 1e-6, silent = FALSE, quiet = TRUE)
}


\details{
 The function \code{SCCcv} is used to perform cross-validation and fit the solution path for the sparse and positive definite basis covariance matrix estimator described in Molstad, Ekvall, and Suder (2024, EJS). The estimator is the minimizer of 
  \deqn{\|\widehat\Theta - \omega 1_p^\top - 1_p \omega^\top + 2 \Omega\|_F^2 + \lambda_1 \sum_{j \neq k}|\Omega_{jk}|}
  subject to \eqn{\Omega \succeq \epsilon I_p}
  where, given \eqn{z_{ijk} = \log(x_{ij}/x_{ik})} and \eqn{\bar{z}_{jk} = \frac{1}{n}\sum_{i=1}^n z_{ijk}}, 
  \deqn{\widehat\Theta_{jk} = \frac{1}{n}\sum_{i=1}^n(z_{ijk} - \bar{z}_{jk})^2.}
  Note that the input data, \code{dat}, are the underlying compositions \eqn{x_{i}= (x_{i1}, \dots, x_{ip})^\top}. For a description of how we select tuning parameters, see the details for \code{gSCCcv}: we use exactly the criterion described therein with \eqn{H=1}. 
}



\arguments{
\item{dat}{Training data: a matrix of dimension \eqn{n \times p}. Each row of the matrix must have positive entries which sum to one.}
\item{nlambda1}{Number of candidate tuning parameters for the lasso-type penalty applied to each basis covariance matrix. }
\item{deltalambda1}{The ratio of minimum to maximum \eqn{\lambda_1} tuning parameter values: must be between zero and one. Decrease if validation errors are minimized for smallest candidate \eqn{\lambda_1} value.}
\item{nfolds}{Positive integer indicating the number of folds to be used for cross-validation.}
\item{alpha0}{Starting step size for proximal gradient descent variants. Leave as default unless you have a compelling reason to modify.}
\item{tau}{Lower bound on eigenvalues of solution. Adjust with caution.}
\item{max.iter}{Maximum number of allowed AccPGD or PPGD iterations. }
\item{tol.obj}{Percent change (divided by 100) in objective function for convergence for AccPGD. Adjust with caution. }
\item{tol.diff}{Maximum absolute difference in successive iterates to claim convergence for AccPGD.Adjust with caution. }
\item{tol.obj.con}{Percent change (divided by 100) in objective function for convergence for PPGD (recommended to be lower than \code{tol.obj}.} Adjust with caution.)
\item{tol.diff.con}{Maximum absolute difference in successive iterates to claim convergence for PPGD. Adjust with caution.}
\item{silent}{TRUE/FALSE; suppress all progress messages? Default is TRUE.}
\item{quiet}{TRUE/FALSE; suppress objective function value printing after each iteration? Should be used only for diagnostic purposes.}
}

\value{
\item{\code{Omega}}{An array of dimension \eqn{1 \times p \times p \times \texttt{nlambda1}}. The first three dimensions correspond to the estimated covariance matrix for the different candidate tuning parameter values.}
\item{\code{valErrs}}{A matrix of dimension \eqn{\texttt{nfolds} \times \texttt{nlambda1}} with validation errors for each candidate tuning parameter.}
\item{\code{lambda1.vec}}{The candidate tuning parameters for \eqn{\lambda_1}. }
\item{\code{avg.valErrs}}{A vector of length \eqn{\texttt{nlambda1}} with the cross-validation errors averaged over the folds.}
\item{\code{Omega.min}}{An array of dimension \eqn{1 \times p \times p} with the covariance matrix estimate using the tuning parameter which minimized the averaged cross-validation error. }
\item{\code{fold.id}}{The fold IDs used for cross-validation.}
}

