% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.fun.R
\name{normalize.max}
\alias{normalize.max}
\title{Max normalization}
\usage{
normalize.max(S)
}
\arguments{
\item{S}{matrix with the raw non normalized scores. Rows are examples and columns are classes}
}
\value{
A score matrix with the same dimensions of \code{S}, but with scores max/normalized separately for each class
}
\description{
Function to normalize the scores of a flat scores matrix per class
}
\details{
The scores of each class are normalized by dividing the score values for the maximum score of that class.
If the max score of a class is zero, no normalization is needed, otherwise \code{NaN} value will be printed as results of 0 out of 0 division.
}
\examples{
data(scores);
maxnorm <- normalize.max(S);
}
