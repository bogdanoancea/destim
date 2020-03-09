#' This function returns the (minus) logLikelihood
#'
#' This function returns the minus logarithm of the likelihood to
#' minimize.
# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first value is expected to have an observation.

logLik <- function(...) {
  UseMethod("logLik")
}
#' @rdname logLik
logLik.HMM <- function(x,e) {
  TM <- createTM(x$transitions, x$parameters$transitions, nstates(x))
  return(floglik(TM, TM@x, rparams(x), x$parameters$reducedparams$transmatrix,
                 istates(x), emissions(x), as.integer(e) - 1L))
}
