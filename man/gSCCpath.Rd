\name{gSCCpath}
\alias{gSCCpath}
\title{gSCCpath}
\description{Fit the solution path for the sparse and positive definite estimator of multiple basis covariance matrices from compositional data}
\usage{
gSCCpath(dat, datval = NULL, nlambda1 = 25, nlambda2 = 25, 
      deltalambda1 = 0.01, deltalambda2 = 0.01, weighted = FALSE,
      lambda1.mat = NULL, lambda2.vec = NULL, alpha0 = 10,  
      tau = 1e-4, max.iter = 5e3, tol.obj = 1e-6, tol.diff = 1e-6,
      tol.obj.con = 1e-6, tol.diff.con = 1e-6, silent = FALSE, quiet = TRUE,
      path.requires = NULL)
}



\arguments{
\item{dat}{Training data: a list of \eqn{H} distinct \eqn{n_{(h)} \times p} matrices for \eqn{h \in \{1, \dots, H\}} with \eqn{H \geq 2}. Each row of the matrices must have positive entries which sum to one.}
\item{datval}{Validation data: a list of \eqn{H} distinct \eqn{n'_{(h)} \times p} matrices for \eqn{h \in \{1, \dots, H\}} with \eqn{H \geq 2}. Each row of the matrices must have positive entries which sum to one.}
\item{nlambda1}{Number of candidate tuning parameters for the lasso-type penalty applied to each basis covariance matrix. }
\item{nlambda2}{Number of candidate tuning parameters for the group-lasso penalty applied to the off-diagonals of the $H$ covariance matrices. }
\item{deltalambda1}{The ratio of minimum to maximum \eqn{\lambda_1} tuning parameter values: must be between zero and one. Decrease if validation errors are minimized for smallest candidate \eqn{\lambda_1} value.}
\item{deltalambda2}{The ratio of minimum to maximum \eqn{\lambda_2} tuning parameter values: must be between zero and one. Decrease if validation errors are minimized for smallest candidate \eqn{\lambda_2} value.}
\item{weighted}{TRUE/FALSE: Should the weighted loss function be used? Default is FALSE.}
\item{alpha0}{Starting step size for proximal gradient descent variants. Leave as default unless you have a compelling reason to modify.}
\item{tau}{Lower bound on eigenvalues of solution. Adjust with caution.}
\item{max.iter}{Maximum number of allowed AccPGD or PPGD iterations. }
\item{tol.obj}{Percent change (divided by 100) in objective function for convergence for AccPGD. Adjust with caution. }
\item{tol.diff}{Maximum absolute difference in successive iterates to claim convergence for AccPGD.Adjust with caution. }
\item{tol.obj.con}{Percent change (divided by 100) in objective function for convergence for PPGD (recommended to be lower than \code{tol.obj}.} Adjust with caution.)
\item{tol.diff.con}{Maximum absolute difference in successive iterates to claim convergence for PPGD. Adjust with caution.}
\item{silent}{TRUE/FALSE; suppress all progress messages? Default is TRUE.}
\item{quiet}{TRUE/FALSE; suppress objective function value printing after each iteration? Should be used only for diagnostic purposes.}
\item{path.requires}{How many tuning parameters must be computed along the solution path? If fewer than \code{nlambda1}, the software will terminate solution path if validation error has increased for five consecutive \eqn{\lambda_1}. NULL (the default) computes the entire path. }
\item{lambda1.mat}{An \code{nlambda2} by \code{nlambda1} matrix of candidate \eqn{\lambda_1}: each column is a vector for candidate \eqn{\lambda_1} for the corresponding value of \eqn{\lambda_2}. NULL allows the software to select tuning parameters internally, so this is recommended. If not NULL, must also specify \code{lambda2.vec}.}
\item{lambda2.vec}{A vector candidate \eqn{\lambda_2} . NULL allows the software to select tuning parameters internally, so this is recommended. }
}

\value{
\item{\code{Omega}:}{An array of dimension \eqn{H \times p \times p \times \texttt{nlambda1} \times \texttt{nlambda2}}. The first three dimensions correspond to the estimated covariance matrices for the \eqn{H} populations for the different candidate tuning parameter values.}
\item{\code{valErrs}:}{A matrix of dimension \eqn{\texttt{nlambda1} \times \texttt{nlambda2}} with validation errors for each candidate tuning parameter pair.}
\item{\code{lambda1.mat}}{The candidate tuning parameters for \eqn{\lambda_1}. Each row corresponds to one \eqn{\lambda_2} value (i.e., for each \eqn{\lambda_2}, we try a slightly different set of \eqn{\lambda_1}). }
\item{\code{lambda2.vec}:}{The candidate tuning parameters for \eqn{\lambda_2}.}
}

