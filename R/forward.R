#' The forward part of the FB algorithm
#'
#' Calculates the forward probabilities.
#'
#' This function has two main purposes: on the one hand it is called by
#' logLik to compute the likelihood, and on the other hand it is combined
#' with backward algorithm to calculate smooth states and smooth
#' consecutive pairwise states by functions sstates and scpstates.
#'
#' @param x A HMM model.
#' @param y A vector with the observed events. It admits missing values.
#'
#' @return A list with two elements:
#' \itemize{
#' \item \emph{alpha} is a matrix that contains the filtered states.
#' The number of columns coincides with the number of observations,
#' so that column k contains the filtered state at the same time as
#' y[k] is observed. For missing (un)observed values, the predicted
#' state is returned instead.
#' \item \emph{scalefactors} is a vector that containts the likelihood
#' of each observation conditioned on all the observation from its
#' past. It is conveniently set to one for missing (un)observed values,
#' so that the joint likelihood is just the cumulative product of the
#' scalefactors.As obvious, its length coincides with the length of
#' y.
#' }
#'
#' @seealso \link{backward}, \link{sstates}, \link{scpstates}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' emissions(model) <- Matrix(c(1, 1, 0.5, 1, 0.5, 0, 0.5, 0, 0,
#'                              0, 0, 0.5, 0, 0.5, 1, 0.5, 1, 1),
#'                              ncol = 2, sparse = TRUE)
#' model <- initparams(model)
#' forward(model, c(1,2,1))

forward <- function(...) {
  UseMethod("forward")
}
#' @rdname forward
forward.HMM <- function(x,y) return(
  fforward(getTM(model), istates(model), emissions(model), as.integer(y) - 1L)
  )

