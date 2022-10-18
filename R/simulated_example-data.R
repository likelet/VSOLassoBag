#' Simulated Example Data for VSOLassoBag Application
#'
#' @docType data
#' @keywords datasets
#' @format ExpressionData is an object constructed by  SummarizedExperiment.
#' it contains the Simulated Example Data for VSOLassoBag with two parts.
#' \describe{
#'     \item{assay(ExpressionData)}{The independent variables matrix (X)
#'     contains 500 variates (rows) x 200 samples (columns).}
#'     \item{colData(ExpressionData)}{The dependent variable(s) matrix (Y)
#'     contains same rows as samples and 1 variate (column) for gaussian,
#'     binomial, possion model application, or 2 variates (columns) for
#'     mgaussian, multinomial and cox model application. The first 1~10
#'     independent variables (X_1~X_10) are simulated to be related to the
#'     dependent variable (D_1), and the first 6~15 independent variables
#'     (X_6~X_15) are simulated to be related to the second dependent variable
#'     (D_2) for mgaussian and multinomial model application. Survival data for
#'     cox model application were simulated with right-censored rate = 0.5
#'     using \code{sim}. \code{survdata} function derived from the coxed R
#'     package.
#'     }
#' }
#'
#' @examples
#' data("ExpressionData")
"ExpressionData"
