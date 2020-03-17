#' Creates a basic rectangle HMM model
#'
HMMrectangle <- function(x,y) {
  x <- as.integer(x)
  y <- as.integer(y)
  if (min(x,y) < 3)
    stop("Both dimensions have to be greater than 2.")
  TL <- createrectangleTL(x,y)
  TL <- TL[, orderTL(TL)]
  CT <- createrectangleCT(TL, x, y)
  newmodel <- HMM(x*y, TL, CT, checks = FALSE)

  return(newmodel)
}
