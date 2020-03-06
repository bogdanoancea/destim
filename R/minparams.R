#' Reparametrizes a HMM model with a minimal set of parameters.
#'
#' Finds a minimal set of parameters that fully determine the
#' probabilities of transition.
#'
#' This function avoids to solve a high dimensional optimization
#' problem with many constraints, parametrizing the probabilities
#' of transition with as few parameters as possible: the number of
#' degrees of freedom.
#'
#' A pivoted QR decomposition of the constraints matrix is done,
#' to get both the free parameters and the matrix that transforms
#' them back into the probabilities of transition.
#'
#' Many constraints are expected to be equalities between two
#' probabilities of transition, so the function is optimized for
#' this special kind of constraints.
#'
#' @param x A HMM object.
#'
#' @return The same HMM object of the input with some additional
#' fields that store the new parameters and the matrix that
#' transforms these parameters in the probabilities of transition.
#'
#' @seealso \link{rparams}, \link{rparams<-}
#'
#' @examples
#' model <- HMMrectangle(2,2)
#' model <- initparams(model)
#' ntransitions(model)
#' nconstraints(model)
#' model <- minparams(model)
#' rparams(model)
#' range(ptransition(model) -
#'   gettransmatrix(model) %*% c(rparams(model), 1))
#'
minparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  reducedparams <- createtransmatrix(x$constraints)
  names(reducedparams)<-c("transmatrix", "numparams")
  if (!is.null(x$parameters$transitions))
    reducedparams$params <- x$parameters$transitions[reducedparams$numparams]
  if (!is.null(x$parameters$reducedparams))
    x$parameters <- x$parameters[c("states", "transitions")]
  x$parameters <- c(x$parameters,list(reducedparams =
                         reducedparams))

  return(x)
  }
