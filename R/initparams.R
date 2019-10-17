#' A simple Initializer for HMM objects
#'
#' Set initial parameters for an HMM object, as specified.
#'
#' The HMM object contains three fields. States contains information
#' about the states. Transitions contains a list of matrices with
#' two rows. Each column represents a transition with non-zero
#' probability. Transitions in the same matrix have same probability.
#' Emissions is a matrix with the emission probabilities. These are
#' considered fixed.
#'
#' @param x A HMM object.
#'
#' @return A initialized HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))
#'
initparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  # Sets random parameters before enforcing constraints
  states <- runif(nstates(x))
  transitions <- runif(ntransitions(x))

  # Enforce constraints: we calculate the closest point in the
  # subspace using Lagrange multipliers
  eqsys <- cbind(matrix(0,nrow = nconstraints(x),
                        ncol = nconstraints(x)),
                 getconstraints(x))
  eqsys <- rbind(eqsys,
                 cbind(t(getconstraints(x)[, 1:ntransitions(x)]),
                       diag(ntransitions(x)) * 2,
                       matrix(2 * transitions, ncol = 1)))
  while (TRUE) {
  transitions <-
    solve(eqsys[,1:nrow(eqsys)], eqsys[,ncol(eqsys)])
  transitions <- transitions[-1:(-nconstraints(x))]
  if (all(transitions >= 0 & transitions <= 1))
    break
  transitions[transitions < 0] <- runif(sum(transitions < 0))
  transitions[transitions > 1] <- runif(sum(transitions > 1))
  eqsys[-1:(-nconstraints(x)), ncol(eqsys)] <- 2 * transitions
  }
  states <- states / sum(states)
  x$parameters <- list(states = states, transitions = transitions)
  return(x)
}
