% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccscore.R
\name{betazeroccscore}
\alias{betazeroccscore}
\title{Hscore beta = 0}
\usage{
betazeroccscore(mat)
}
\arguments{
\item{mat}{a matrix}
}
\value{
The hscore of \code{mat}
}
\description{
Calculates Hscore with beta equal to 0
}
\examples{
mat <- replicate(5, rnorm(5))
betazeroccscore(mat)
}
