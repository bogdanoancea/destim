#' @title Initializer for HMM objects
#'
#' @description Sets initial parameters for a HMM object, as specified.
#' The field parameters of the HMM object, which includes both
#' initial state and transition probabilities, is initialized at
#' random.
#'
#' The initial states probabilities are set to an uniform
#' (0,1) distribution and then divided by their sum.
#'
#' The initial probabilities of transition are also set to an uniform
#' (0,1) and in this case, projected on the constrained space. After
#' the projection some probability might result greater than one or
#' less than zero. Those probabilities are then set to uniform (0,1)
#' again and the process is repeated until all probabilities of
#' transition are in (0,1) and the constraints are satisfied.
#'
#' @param x A HMM object.
#'
#' @return An initialized HMM object.
#'
#' @seealso \link{HMM}, \link{minparams}, \link{initsteady}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' model <- initparams(model)
#' range(constraints(model) %*% c(ptransition(model), -1)) # It should be close to zero
#'
#' @export
initparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  # Sets initial state at random
  states <- runif(nstates(x))
  states <- states / sum(states)

  x$parameters <- list(states = states,
                       transitions = inittransitions(x$constraints))
  return(x)
}
