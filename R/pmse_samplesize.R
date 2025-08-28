#' pmse_samplesize
#' - Sample Size Calculation for Prediction Models
#'
#' @title Compute efficient sample size under user-defined PMSE targets
#'
#' @description \code{pmse_samplesize} computes a sample size for a
#' prediction model. The function implements the formulas found in the thesis
#' "Predictive Power and Efficient Sample Size in Linear Regression Models" by Yifan Ma (2023).
#'
#' @details \code{pmse_samplesize} The function calculates predictor error variance
#' for the full model, with all predictors, and the reduced model, with the basic
#' predictors using a provided covariance matrix or correlation matrix. It can
#'  also calculate predictor error variance through Cohen's F^2 and R^2 values.
#'  With the predictor error variance it determines a sample size from the
#'  efficient sample size at a target efficiency level and a sample size from a
#'  PMSE value of the full and reduced model. The final returned sample size is
#'  the largest out of the outputs.
#'
#'
#' @param k Integer. Total number of predictors in the full model.
#'
#' @param p Integer. Number of basic predictors in the reduced model.
#'
#' @param PMSE_val_k Numeric. Target PMSE value for the full model.
#'
#' @param PMSE_val_p Numeric. Target PMSE value for the reduced model.
#'
#' @param efficiency_level Numeric. Target efficiency level.
#' (default is 0.9, meaning 90% of asymptotic pPMSEr)
#'
#' @param sigma_k2 Numeric. Predictor error variance for full model. If 'NULL'
#' it is derived.
#'
#' @param sigma_p2 Numeric. Predictor error variance for basic model. If 'NULL'
#' it is derived.
#'
#' @param cov Optional covariance matrix. Must be `(k+1) x (k+1)` with the response
#' 1st row and column.
#'
#' @param corr Optional correlation matrix. (Same layout as `cov`).
#'
#' @param SD Optional numeric vector of standard deviation for the predictors when
#' a correlation matrix is supplied. Default `1`
#'
#' @param f2 Numeric. Cohen's f2 for effects of all predictors in full model.
#'
#' @param f2_2 Numeric. Cohen’s f2 for the effects of new predictors given
#' the basic model.
#'
#' @param R2_full Numeric. Coefficient of determination for full model.
#'
#' @param R2_basic Numeric. Coefficient of determination for basic model.
#'
#' @return Numeric representing the required sample size.
#'
#' @references
#' Ma, Y. (2023). _Predictive Power and Efficient Sample Size in Linear
#' Regression Models_. Master’s Thesis, Worcester Polytechnic Institute.
#'
#' @examples
#' ## Example with a 5-predictor model (k = 5) and 2 basic predictors (p = 2)
#' pmse_samplesize(
#'   k = 5, p = 2,
#'   PMSE_val_k    = 1,
#'   PMSE_val_p    = 1,
#'   efficiency_level = 0.9,
#'   sigma_k2 = 0.50,
#'   sigma_p2 = 0.60
#' )
#'
#' @export
pmse_samplesize = function(k, p,
                           PMSE_val_k = 1,
                           PMSE_val_p = 1,
                           efficiency_level = 0.9,
                           sigma_k2 = NULL,
                           sigma_p2 = NULL,
                           cov = NULL,
                           corr = NULL,
                           SD = 1,
                           f2 = NULL,
                           f2_2 = NULL,
                           R2_full = NULL,
                           R2_basic = NULL) {

  sample_size <- -1
  #Make sure p and k are entered
  # If there is a covariance matrix, check if it has the right dimensions


  # Get covariance matrix if isn't given
  if(is.null(cov) & !is.null(corr)){
    cov <- Corr_to_Cov(corr, SD)
  }

  # Get sigma k and p values with covariance matrix
  if(!is.null(cov) & is.null(sigma_k2)){
    sigma_k2 <- sigmak2(cov)
  }
  if(!is.null(cov) & is.null(sigma_p2)){
    sigma_p2 <- sigmap2(cov, p, k)
  }

  # Get sigma k and p from R^2 and F^2 values
  # Obtain F^2 from R^2 if it does not exist
  if(is.null(f2_2) & !is.null(R2_full) & !is.null(R2_basic)){
    f2_2 <- f2_2_calcu(R2_full, R2_basic)
  }
  # Returns only efficient sample size if only given F2 and R2 values.
  if(!is.null(f2_2) & is.null(sigma_k2) & is.null(sigma_p2)){
    sigma_p2 <- f2_2 + 1
    simga_k2 <- 1
    effN <- gen_efficient_sampsize(1, sigma_p2, k, p, 1 - efficiency_level)
    return(effN)
  }

  # with sigma_k and sigma_p values, obtain sample size from PMSE and efficient sample size
  if(!is.null(sigma_k2) & !is.null(sigma_p2)){
    pmseK <- gen_PMSE_sampsize(sigma_k2, k, PMSE = PMSE_val_k)
    pmseP <- gen_PMSE_sampsize(sigma_p2, p, PMSE = PMSE_val_p)
    effN <- gen_efficient_sampsize(sigma_k2, sigma_p2, k, p, 1 - efficiency_level)
    if(pmseK == -1){
      stop("PMSE value for full model is too small")
    }
    if(pmseP == -1){
      stop("PMSE value for basic model is too small")
    }
    sample_size <- max(sample_size, pmseK, pmseP, effN)
    return (sample_size)
  }
  else{
    stop("Can not determine predictor error variance")
  }
}





