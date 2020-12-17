#' @title Minus logLikelihood
#'
#' @description Returns the minus logarithm of the likelihood given a model and a set of observations.
#'
#' A slightly modified version of the forward algorithm is used to compute the likelihood,
#' to avoid store unneeded data. The sign is changed because it is usual to minimize
#' instead maximize.
#'
#' @param x A HMM model.
#' @param e A vector with the observed events. It admits missing values.
#'
#' @return The minus logarithm of the likelihood of the events given the model.
#'
#' @seealso \link{fit}, \link{HMM}, \link{initparams}
#'
#' @examples
#' model <- HMMrectangle(20,20)
#' S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
#' emissions(model) <- createEM(c(20,20), towers, S)
#' model <- initparams(model)
#' model <- minparams(model)
#' logLik(model,events)
#'
#' @export
logLik <- function(...) {
  UseMethod("logLik")
}
#' @rdname logLik
#' @export
logLik.HMM <- function(x,e, ...) {
  TM <- createTM(x$transitions, x$parameters$transitions, nstates(x))
  return(floglik(TM, TM@x, rparams(x), x$parameters$reducedparams$transmatrix,
                 istates(x), emissions(x), as.integer(e) - 1L))
}
