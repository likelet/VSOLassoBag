% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpleEstimation.R
\name{simpleEstimation}
\alias{simpleEstimation}
\title{Parametric Statistical Test}
\usage{
simpleEstimation(res.df, bootN)
}
\arguments{
\item{res.df}{a dataframe with variables and observed frequency}

\item{bootN}{an integer, bagging times}
}
\value{
a list of p-value of each variable and the average selection ratio
}
\description{
An internal function utilized by VSOLassoBag.
}
\examples{
load(system.file("extdata/Results.RData", package="VSOLassoBag"))
simpleEstimation(Results, 10)
}
\references{
\href{https://www.doi.org/10.1371/journal.pone.0141869}{RRLASSO, Park H., et al, 2015}, the algorithm utilized in LassoBag has been modified.
}
