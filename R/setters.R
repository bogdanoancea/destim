#' Functions to modify a HMM object
#'
`rparams<-` <- function(x, value) {
  if (is.null(x$parameters$reducedparams)) {
    if (is.null(x$parameters$states))
      x <- initparams(x)
    x <- minparams(x)
  }
  x$parameters$reducedparams$params <- value
  x$parameters$transitions <-
    (x$parameters$reducedparams$transmatrix %*%
      matrix(c(x$parameters$reducedparams$params, 1), ncol = 1)) [,]
  return(x)
}
