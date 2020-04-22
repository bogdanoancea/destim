#' The backward part of the FB algorithm
#'
#' Calculates the backward probabilities.
#'
#' The main purpose of this function is to be combined with forward function
#' to calculate smooth states and smooth consecutive pairwise states. This is done
#' by functions sstates and scpstates.
#'
#' @param x A HMM model.
#' @param y A vector with the observed events. It admits missing values.
#' @param sfactors A vector with the scale factors from the forward algorithm.
#'
#' @return A matrix that contains the (rescaled) backward probabilities.
#'
#' @seealso \link{forward}, \link{sstates}, \link{scpstates}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' emissions(model) <- Matrix(c(1, 1, 0.5, 1, 0.5, 0, 0.5, 0, 0,
#'                              0, 0, 0.5, 0, 0.5, 1, 0.5, 1, 1),
#'                              ncol = 2, sparse = TRUE)
#' model <- initparams(model)
#' fpass <- forward(model, c(1,2,1))
#' backward(model, c(1,2,1), fpass$scalefactors)

backward <- function(...) {
  UseMethod("backward")
}
#' @rdname backward
backward.HMM <- function(x,y,sfactors) return (
  fbackward(getTM(x), emissions(x), as.integer(y) - 1L, sfactors)
  )
