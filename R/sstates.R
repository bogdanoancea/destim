#' Smooth states
#'
#' Returns the smooth states from the forward-backward algorithm.
#'
#' Smooth states are the marginal of the state conditional on the observations, for each time.
#' This agrees with the so called \eqn{\gamma} from the Baum-Welch algorithm.
#'
#' It is returned as a matrix, so that the smooth state for time instant i is the column i of
#' the matrix.
#'
#' @param x A HMM model.
#' @param e A vector with the observed events. It admits missing values.
#'
#' @return A sparse matrix. The number of rows is the number of states, and the number of columns
#' is the number of observed events. Each column of the output matrix corresponds to the
#' probability mass function for the state, so it sums up to one.
#'
#' @seealso \link{HMM}, \link{scpstates}, \link{backward}
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
#' sstates(model, obs)
#'
#' @export
sstates <- function(...) {
  UseMethod("sstates")
}
#' @rdname sstates
#' @export
sstates.HMM <- function(x, e) {
  fpass <- forward(x,e)
  bpass <- backward(x,e,fpass$scalefactors)
  return(fpass$alpha * bpass)
}
