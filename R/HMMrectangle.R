#' Creates a basic rectangle HMM model
#'
HMMrectangle <- function(x,y) {
  x <- as.integer(x)
  y <- as.integer(y)
  if ((x <= 0) || (y <= 0))
    stop("Both dimensions have to be integer positive numbers.")
  TL <- createrectangleTL(x,y)
  CT <- createrectangleCT(TL, x, y)
  newmodel <- HMM(x*y, TL, CT, checks = FALSE)

  return(newmodel)
}
