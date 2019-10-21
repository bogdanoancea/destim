#' Class constructor for hidden Markov models
#'
#' Creates an HMM object, as specified.
#'
#' The HMM object contains three fields. States contains information
#' about the states. Transitions contains a list of matrices with
#' two rows. Each column represents a transition with non-zero
#' probability. Transitions in the same matrix have same probability.
#' Emissions is a matrix with the emission probabilities. These are
#' considered fixed.
#'
#' @param S  Number or list of states.
#' @param TL Matrix with non-zero transitions.
#' @param CT Matrix of constraints.
#' @param EM Matrix of emissions.
#'
#' @return A HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))

HMM <- function(...) {
  UseMethod("HMM")
}
#' @rdname HMM
HMM.integer <- function(S, TL, CT, EM = NULL) {
  if (missing(TL))
    TL = matrix(c(1:S,1:S), nrow = 2, byrow = TRUE)
  if (missing(CT))
    CT <- t(sapply(1:S, function(x) c(TL[1,] == x, 1.0)))
  output <- list(states = list(names = as.character(1:S),
                               coordinates = NULL),
                 transitions = TL, constraints = CT,
                 emissions = EM, parameters = list(transitions = NULL,
                                                   states = NULL))
  attr(output,"class") <- "HMM"
  return(output)
}
#' @rdname HMM
HMM.character <- function(S, ...) {
  output <- HMM(length(S), ...)
  setsnames(output) <- S
  return(output)
}
