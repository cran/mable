#' Pancreatic Cancer Biomarker Data
#'
#' @description Contain sera measurements from 51 control patients with pancreatitis and 90
#' case patients with pancreatic cancer at the Mayo Clinic with a cancer antigen, CA125,
#'  and with a carbohydrate antigen, CA19-9 (Wieand, et al, 1989)
#' @source Wieand, S., Gail, M. H., James, B. R., and James, K.L. (1989). 
#'  A family of nonparametric statistics for comparing diagnostic markers with 
#'  paired or unpaired data. Biometrika, 76, 585--592.
#'
#' @docType data
#' @keywords datasets
#' @name pancreas
#' @usage data(pancreas)
#' @format A data frame with 141 rows and 3 variables. 
#' 
#' \itemize{
#'   \item ca199. CA19-9 levels
#'   \item ca125. CA125 levels
#'   \item status. 0 = controls (non-cancer) and 1 = cases (cancer).
#' }
#' @references 
#' Wieand, S., Gail, M. H., James, B. R., and James, K.L. (1989). 
#'  A family of nonparametric statistics for comparing diagnostic markers with 
#'  paired or unpaired data. Biometrika, 76, 585--592.
#'
#' @examples data(pancreas)
"pancreas"
