#' Convert a correlation matrix to a covariance matrix
#'
#' Computes \eqn{\Sigma = D\,\mathrm{CORR}\,D} where
#' \eqn{D = \mathrm{diag}(\texttt{SD})}.
#'
#' @param CORR Square correlation matrix.
#' @param SD   Numeric vector of standard deviations (same length as `ncol(CORR)`).
#'
#' @return Covariance matrix of identical dimension to `CORR`.
#' @examples
#' C <- matrix(c(1, .3, .3, 1), 2)
#' Corr_to_Cov(C, SD = c(2, 1))
#' @keywords internal
#' @noRd
Corr_to_Cov <- function(CORR, SD) {
  sweep(sweep(CORR, 1, SD, "*"), 2, SD, "*")
}

#' Error-term variance for the full regression model
#'
#' Calculates \eqn{\sigma_k^2 = \sigma_{00} -
#'   \boldsymbol{\sigma}^\top \Sigma^{-1}\boldsymbol{\sigma}}.
#'
#' @param COV Covariance matrix whose first row/column correspond to the
#'   response followed by the \eqn{k} predictors.
#'
#' @return Numeric scalar (`sigma_k2`).
#' @keywords internal
#' @noRd
sigmak2 <- function(COV){
  sigma_00 <- COV[1, 1]
  sig      <- as.matrix(COV[-1, 1])
  SIG      <- COV[-1, -1]
  sigma_00 - t(sig) %*% solve(SIG) %*% sig
}

#' Error-term variance for the reduced regression model
#'
#' Computes \eqn{\sigma_p^2} when only the first \eqn{p} of \eqn{k}
#' predictors are retained.
#'
#' @param COV Covariance matrix of the response and all \eqn{k} predictors.
#' @param p   Number of predictors in the reduced model.
#' @param k   Total number of predictors in the full model (unused inside).
#'
#' @return Numeric scalar (`sigma_p2`).
#' @keywords internal
#' @noRd
sigmap2 <- function(COV, p, k){
  sigma_00 <- COV[1, 1]
  sig      <- as.matrix(COV[-1, 1])
  SIG      <- COV[-1, -1]
  sig_1    <- sig[1:p, , drop = FALSE]
  SIG_11   <- SIG[1:p, 1:p, drop = FALSE]
  sigma_00 - t(sig_1) %*% solve(SIG_11) %*% sig_1
}

#' Cohen’s *f*² from two *R*² values
#'
#' Computes \eqn{f^2 = (R^2 - R_1^2)/(1 - R^2)}.
#'
#' @param R2   *R*-squared for the full model.
#' @param R2_1 *R*-squared for the reduced model.
#'
#' @return Numeric scalar (`f2`).
#' @keywords internal
#' @noRd
f2_2_calcu <- function(R2, R2_1){
  (R2 - R2_1) / (1 - R2)
}

#' Efficient sample size for a target pPMSEr fraction
#'
#' Finds the minimum \eqn{n} that achieves a fraction `alpha`
#' of the maximal predictive-MSE reduction.
#'
#' @param sigma_k2 Error variance of the full model.
#' @param sigma_p2 Error variance of the reduced model.
#' @param k        Number of predictors in the full model.
#' @param p        Number of predictors in the reduced model.
#' @param alpha    Fraction of the maximum pPMSEr.
#'
#' @return Numeric scalar (`n_star`).
#' @keywords internal
#' @noRd
gen_efficient_sampsize <- function(sigma_k2, sigma_p2, k, p, alpha = 0.1){
  EVR         <- sigma_k2 / sigma_p2
  if (!is.finite(EVR) || EVR >= 1) {
    stop("No efficiency gain (sigma_k2 >= sigma_p2).")
  }
  lambda_star <- 1 + alpha * (1 / EVR - 1)
  p2          <- k - p
  p + 2 + (p2 * lambda_star / (lambda_star - 1))
}

#' Sample size required to achieve a target PMSE
#'
#' Uses `rootSolve::uniroot.all()` to find the smallest \eqn{n}
#' such that the predictive mean-square error equals `PMSE`.
#'
#' @param variance Error variance (\eqn{\sigma^2}).
#' @param k        Number of predictors in the model.
#' @param PMSE     Target PMSE.
#'
#' @return Numeric vector of feasible sample sizes,
#'   or \code{-1} if the target is unattainable (`PMSE <= variance`).
#' @importFrom rootSolve uniroot.all
#' @keywords internal
#' @noRd
gen_PMSE_sampsize <- function(variance, k, PMSE){
  variance   <- as.vector(variance)
  PMSE_fun <- function(n){
    variance * (n + 1) * (n - 2) / (n * (n - k - 2)) - PMSE
  }
  if (PMSE <= variance) return(-1)

  # Exponential search for an upper bound
  ub <- 300
  while (PMSE_fun(ub) >= 0) {
    ub <- ub * 2
    if (ub > 1e6)
      stop("Upper bound too high; check inputs.")
  }
  rootSolve::uniroot.all(PMSE_fun, c(k + 2, ub))
}

