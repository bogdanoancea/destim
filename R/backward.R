#' The backward part of the FB algorithm
#'
#' Calculates the backward probabilities.
#'
#'
#' @return A HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))

backward <- function(...) {
  UseMethod("backward")
}
#' @rdname backward
backward.HMM <- function(x,y,sfactors) return (
  fbackward(getTM(model), emissions(model), as.integer(y) - 1L, sfactors)
  )
