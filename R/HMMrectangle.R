#' Basic HMM grid model
#'
#' Creates a basic rectangular grid model as specified.
#'
#' This model is a rectangular grid where the only transitions allowed are those between
#' contiguous tiles. Moreover, horizontal and vertical transition probabilities are
#' equal for all tiles. Diagonal transition probabilities are equal between them too,
#' but different from the former. These constraints mean that there are only two parameters
#' to estimate.
#'
#' The emissions field is left unasigned.
#'
#' @param x length of the rectangle in tiles.
#' @param y width of the rectangle in tiles.
#'
#' @return A HMM object.
#'
#' @seealso \link{emissions}, \link{minparams}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' nstates(model)
#' ntransitions(model)
#' nconstraints(model)
HMMrectangle <- function(x,y) {
  x <- as.integer(x[1])
  y <- as.integer(y[1])
  if (min(x,y) < 3)
    stop("Both dimensions have to be greater than 2.")
  TL <- createrectangleTL(x,y)
  TL <- TL[, orderTL(TL)]
  CT <- createrectangleCT(TL, x, y)
  newmodel <- HMM(x*y, TL, CT, checks = FALSE)

  return(newmodel)
}
