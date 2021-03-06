% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setters.R
\name{emissions<-}
\alias{emissions<-}
\title{Set the emissions of the model}
\usage{
emissions(x) <- value
}
\arguments{
\item{x}{A HMM model.}

\item{value}{A (sparse column major) matrix with the likelihoods of each emission (column)
conditioned on the state (row).}
}
\value{
Changes the emissions matrix in the model.
}
\description{
Sets the emissions of the model.

The number of rows must match the number of states. If a matrix is provided,
it is converted to column major sparse matrix (\code{dgCMatrix}).
}
\examples{
model <- HMMrectangle(10,10)
tws <- matrix(c(3.2, 6.1, 2.2, 5.7, 5.9, 9.3, 5.4,
4.0, 2.9, 8.6, 6.9, 6.2, 9.7, 1.3),
nrow = 2, ncol = 7)
S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
emissions(model)<-createEM(c(10,10), tws, S)
emissions(model)

}
\seealso{
\link{HMM}, \link{emissions}
}
