#' @title Creates the events matrix.
#'
#' @description   Creates the events matrix for a rectangular grid
#' according to the location of the towers and the S function.
#'
#' @param size The number of rows and columns in the grid.
#'
#' @param towers The towers(antenna) positions. This parameter is a matrix with
#' two rows and the number of columns equals to the number of antennas. On the
#' first row we have the X coordinate of the towers and on the secnd row the
#' Y coordinate.
#'
#' @param S A function to compute the signal strength.
#'
#' @return The events matrix.

#' @export
createEM <- function(size, towers, S) {
  output <- matrix(0, nrow = size[1] * size[2],
                   ncol = ncol(towers))

  for (i in 1:size[1])
    for (j in 1:size[2])
      output[(j - 1) * size[2] + i, ] <-
        getER(c(i,j), towers, S)

  return(output)
}
