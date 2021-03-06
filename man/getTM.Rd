% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters.R
\name{getTM}
\alias{getTM}
\alias{getTM.HMM}
\title{Transition matrix.}
\usage{
getTM(x)

\method{getTM}{HMM}(x)
}
\arguments{
\item{x}{the HMM object.}
}
\value{
A row major sparse matrix which is the transition matrix of the model.
}
\description{
Returns the transition matrix from a HMM object.
The transition matrix is represented as row major. This way its transpose matrix
which is used to left multiply the state is column major.
}
\examples{
model <- HMM(2)
model <- addtransition(model,c(1,2))
model <- initparams(model)
getTM(model)

}
\seealso{
\link{HMM}, link{emissions}
}
