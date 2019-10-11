#' Sets the initial state to the steady state
#'
#' The initial a-priori distribution is set to the steady state of
#' the transition matrix.
#'
#' @param x A HMM object.
#'
#' @return A initialized HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))
#'
initsteady <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  # Initial steady state calculation
  edecomp <- eigen(getTM(x))
  # The steady state might be not unique! Fix this later
  i <- which(abs(edecomp$values - 1) < 1e-6)[1]
  states <- Re(edecomp$vectors[,i])
  states <- states / sum(states)

  x$parameters$states <- states

  return(x)
}
