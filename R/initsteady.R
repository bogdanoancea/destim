#' Sets the initial state to the steady state
#'
#' The initial a priori distribution is set to the steady state of
#' the transition matrix.
#'
#' The Markov Chain is expected to be irreducible and aperiodic. The
#' first because otherwise the devices would not have freedom of
#' movement. The second because some probabilities from one state
#' to itself are expected to be non zero. This implies that there
#' exists one unique steady state.
#'
#' The steady state is computed by solving the sparse linear system (TM - I)x = 0, where
#' TM is the matrix of transitions I is identity and x the steady state. As it is an
#' homogeneous system, and because of the uniqueness of the steady state, the solution is
#' a one dimensional vector space, and the generator does not have any coordinate equal to zero.
#' Then the last coordinate is set to 1 / number of states, so the sparse linear system
#' becomes inhomogeneous with unique solution. Finally the solution is normalized so that
#' the components of x sum up to 1.
#'
#' @param x A HMM object.
#'
#' @return The same HMM object of the input with its initial state set
#' to steady state.
#'
#' @seealso \link{HMM}, \link{initparams}, \link{minparams}
#'
#' @examples
#' model <- HMM(2)
#' model <- addtransition(model, c(1,2))
#' model <- addtransition(model, c(2,1))
#' model <- initparams(model)
#' istates(model)
#' model <- initsteady(model)
#' istates(model)
#' (istates(model) %*% getTM(model))
#'
#' @export
initsteady <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  istates(x) <- createsteady(getTM(x))
  return(x)
}
