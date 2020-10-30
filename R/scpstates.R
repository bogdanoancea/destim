#' @title Returns \eqn{\xi} like in the Baum-Welch algorithm.
#'
#' @description Returns the smooth joint probability mass function for consecutive
#' states, which is usually called \eqn{\xi} in the Baum-Welch algorithm.
#' Smooth states are marginal but as they are far to be independent it is convenient to
#' have some information about their dependence. This function returns the joint probability
#' mass function for two time consecutive states, conditional on the observations. This agrees
#' with the so called \eqn{\xi} from the Baum-Welch algorithm.
#'
#' It is returned as a matrix, so that the said joint probability for time instants i - 1 and i
#' are the columns from i - 1 times the number of states plus one, to i times the number of
#' states.
#'
#' @param x A HMM model.
#' @param e A vector with the observed events. It admits missing values.
#'
#' @return A sparse matrix. The number of rows is the number of states, and the number of columns
#' is the number of states times the number of observed events minus one. Each full row square
#' slice of the output matrix corresponds to a joint probability mass function, so it sums up
#' to one.
#'
#' @seealso \link{HMM}, \link{sstates}, \link{backward}
#'
#' @examples
#' model <- HMMrectangle(10,10)
#' tws <- matrix(c(3.2, 6.1, 2.2, 5.7, 5.9, 9.3, 5.4,
#' 4.0, 2.9, 8.6, 6.9, 6.2, 9.7, 1.3),
#' nrow = 2, ncol = 7)
#' S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
#' emissions(model)<-createEM(c(10,10), tws, S)
#' obs <- c(1,2,NA,NA,NA,NA,7,7)
#' model <- fit(model, obs)
#' scpstates(model, obs)
#' @export
scpstates <- function(...) {
  UseMethod("scpstates")
}
#' @rdname scpstates
#' @export
scpstates.HMM <- function(x, e) {
  TM <- getTM(x)
  EM <- emissions(x)
  fpass <- forward(x,e)
  bpass <- backward(x,e,fpass$scalefactors)

  return(
    fscpstates(TM, fpass$alpha, bpass, fpass$scalefactors, EM, as.integer(e) - 1L)
  )
}
