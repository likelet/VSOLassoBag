# provide p-value using the statistical method described in RRLASSO (Park H., et al, 2015; doi:10.1371/journal.pone.0141869)

# Return a p-value list and the average selection ratio pi

#' Parametric Statistical Test
#'
#' This is an internal function utilized by VSOLassoBag.
#'
#' @param res.df a dataframe with variables and observed frequency
#' @param bootN an integer, bagging times
#' @return a list of p-value of each variable and the average selection ratio
#' @references \href{https://www.doi.org/10.1371/journal.pone.0141869}{RRLASSO, Park H., et al, 2015}, the algorithm utilized in LassoBag has been modified.
#' @export
#'
#' @examples
#' data("simulated_example")
#' simpleEstimation(test.df$result,10)

utils::globalVariables("pbinom")
simpleEstimation<-function(res.df,bootN){
  pin<-sum(res.df$Frequency)/(bootN*nrow(res.df))
  pvalue.list<-pbinom(q=res.df$Frequency,size=bootN,prob=pin,lower.tail=FALSE)
  return(list(pvalue=pvalue.list,pi=pin))
}
