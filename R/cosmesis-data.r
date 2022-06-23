#' Breast cosmesis data
#'
#' @description Data contain the interval-censored times to cosmetic deterioration  
#'     for breast cancer patients undergoing radiation or radiation plus chemotherapy.
#' @source Finkelstein, D. M. and Wolfe, R. A. (1985) 
#'     A semiparametric model for regression analysis of interval-censored failure time data. 
#'     Biometrics 41, 933–945.
#'
#' @docType data
#' @keywords datasets
#' @name cosmesis
#' @usage data(cosmesis)
#' @format A data frame with 94 observations on the following 3 variables. 
#' 
#' \itemize{
#'   \item left left endpoint of the censoring interval in months
#'   \item right right endpoint of the censoring interval in months
#'   \item treat a factor with levels \code{RT} and \code{RCT} representing radiotherapy-only 
#'          and radiation plus chemotherapy treatments, respectively
#' }
#' @references 
#' Finkelstein, D. M. (1986) A proportional hazards model for interval-censored 
#'    failure time data. Biometrics 42, 845–854.
#'
#' @examples data(cosmesis)
"cosmesis"
