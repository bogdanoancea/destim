#' Adds constraints to the model
#'
#' The specified contraints are added to the model.
#'
#' If parameter ct is a vector, it is expected to be a set of transition
#' probabilities indexed as in field transitions of the model. In this
#' case the constraint added is the equality between the refered probabilities of
#' transition.
#'
#' If parameter ct is a matrix, it is expected to be a system of additional
#' linear equalities that the model must fulfill. Thus, the new equations
#' are added to the field constraints of the model.
#'
#' While it is possible to use a matrix to add equality constraints, it is not recommended
#' because of performance.
#'
#' Previous constraints of the model are preserved.
#'
#' @param x A HMM object.
#' @param ct The additional constraints, which can be either a
#' matrix or a vector (see details).
#'
#' @return A HMM object similar to the input but with the additional
#' constraints.
#'
#' @seealso \link{HMM}, \link{addtransition}
#'
#' @examples
#' model <- HMM(3)
#' model <- addtransition(model, c(1,2))
#' model <- addtransition(model, c(2,3))
#' model <- addtransition(model, c(3,1))
#' transitions(model)
#' constraints(model)
#' model <- addconstraint(model,c(2,4,5))
#' constraints(model)
#'
#' @export
addconstraint <- function(x, ct) {

  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  if (is.vector(ct)) {
    if (any(ct < 1) || any(ct > ntransitions(x)))
      stop(paste0("The indexes must be numbers between one and the number ",
                  "of transitions."))
    ct <- as.integer(ct)
    stillt <- which(transitions(x)[1,] == transitions(x)[2,]) - 1L
    newct <- faddconstraint(ct - 1L, transitions(x), constraints(x), stillt, nstates(x))

    x[["constraints"]] <- newct
  }
  else if (is.matrix(ct) || (class(ct) == "dgCMatrix")) {
    if (ncol(ct) != (ntransitions(x) + 1))
      stop(paste0("The number of columns of the constraints must ",
                  "agree with the number of transitions plus one."))
    if (is.matrix(ct))
      ct <- Matrix::Matrix(ct, sparse = TRUE)
    x[["constraints"]] <- rbind(x[["constraints"]], ct)
  }

  x[["parameters"]] <- list(transitions = NULL, states = NULL)

  return(x)
}
