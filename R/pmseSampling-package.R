#' pmsesampling: Sample Size Determination for Accurate Predictive Linear Regression
#'
#' Tools to estimate the minimum sample size required to achieve a
#' target Prediction Mean-Squared Error (PMSE) or a specified
#' proportional PMSE reduction (pPMSE<sub>r</sub>).  Functions implement the
#' analytic and simulation-based criteria described in
#' Ma (2023) and include helpers for covariance-matrix handling,
#' root-finding and diagnostic plotting.
#'
#' @section Core functions:
#' \describe{
#'   \item{\code{pmse_samplesize()}}{Determines sample size from PMSE equation in
#'   basic and full models and the efficient sample size}
#' }
#'
#' @section Typical workflow:
#' 1. Obtain \eqn{\sigma_k^2} and \eqn{\sigma_p^2}
#' 2. Or import or build a predictor covariance matrix.\cr
#' 3. Or obtain \eqn{Cohen's f^2} and \eqn{\R^2}
#' 4. Call \code{pmse_samplesize} with available inputs to get sample size.\cr
#' @references
#' Ma Y. (2023) *Predictive Power and Efficient Sample Size in
#' Linear Regression Models*. Worchester Polytechnic Institute
#'
#' @docType package
#' @name pmsesampling-package
#'
#' @importFrom stats uniroot rnorm lm
#' @importFrom Matrix forceSymmetric
"_PACKAGE"
