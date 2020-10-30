#' @title This function returns the viterbi path.
#'
#' @description This function implements the viterbi algorithm to return the
#' most likely path.
#'
#' @param x The HMM model.
#' @param e Sequence of observation events. Since we are taking
#' time increase small, it is expected to have mostly missing
#' values. The first value is expected to have an observation.
#' @return The viterbi path.
#' @keywords internal
viterbi <- function(...) {
  UseMethod("viterbi")
}
#' @rdname viterbi
viterbi.HMM <- function(x,e) {
  output <- integer(length(e))
  TM <- getTM(x)
  EM <- emissions(x)
  W <- matrix(0, nrow = nstates(x), ncol = length(e))
  WI <- matrix(0L, nrow = nstates(x), ncol = length(e))
  W[, 1] <- log(istates(x))
  if (!is.na(e[1]))
    W[, 1] <- W[, 1] + log(EM[,e[1]])
  for (i in 2:length(e)) {
    WI[, i] <- apply(W[, i - 1] + log(TM), 2, which.max)
    W[, i] <- apply(W[, i - 1] + log(TM), 2, max)
    if (!is.na(e[i]))
      W[, i] <- W[, i] + log(EM[,e[i]])
  }
  output[length(e)] <- which.max(W[, length(e)])
  print(max(W[,length(e)]))
  for (i in length(e):2) {
    output[i - 1] <- WI[output[i], i]
  }
  return(output)
}
