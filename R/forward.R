#' The forward part of the FB algorithm
#'
#' Calculates the forward probabilities.
#'
#'
#' @return A HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))

forward <- function(...) {
  UseMethod("forward")
}
#' @rdname forward
forward.HMM <- function(x,y) {
  TM <- getTM(x)
  EM <- emissions(x)
  alpha <- matrix(0,nrow = nstates(x), ncol = length(y))
  sfactors <- numeric(length(y))
  svector <- getiparameters(x)
  for (i in 1:length(y)) {
    if (!is.na(y[i])) {
      svector <- svector * EM[,y[i]]
      sfactors[i] <- sum(svector)
      svector <- svector / sum(svector)
    }
    else sfactors[i] <- 1
    alpha[,i] <- svector
    svector <- svector %*% TM
  }
  return(list(alpha = alpha, scalefactors = sfactors))
}
