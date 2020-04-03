#' Adds constraints to the model
#'
#' The specified contraints are added to the model.
#'
#' If parameter ct is a matrix, it is expected to be a system of additional
#' linear equalities that the model must fulfill. Thus, the new equations
#' are added to the field constraints of the model.
#' If parameter ct is a vector, it is expected to be a set of transition
#' probabilities indexed as in field transitions of the model. In this
#' case the constraint added is the equality between the refered probabilities of
#' transition.
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
addconstraint <- function(x, ct) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  if (is.matrix(ct)) {
    if (ncol(ct) != (ntransitions(x) + 1))
      stop(paste0("The number of columns of the constraints must ",
           "agree with the number of transitions plus one."))
    x[["constraints"]] <- rbind(x[["constraints"]], ct)
  }
  else {
    if (any(ct < 1) || any(ct > ntransitions(x)))
      stop(paste0("The indexes must be numbers between one and the number ",
                  "of transitions."))
    ct <- as.integer(ct)
    print("Paso 1")
    bct <- createBCT(transitions(x), nstates(x))
    print("Paso 2")
    newct <- createEQ(ct - 1L, ntransitions(x) + 1L)
    stillt <- which(transitions(x)[1,] == transitions(x)[2,]) - 1L
    print("Paso 3")
    oldeq <- extractEQ(constraints(x), stillt)
    print("Paso 4")
    if (nrow(oldeq) != 0L)
      newct <- frbind(oldeq, newct)
    print("Paso 5")
    newct <- canonEQ(newct)
    print("Paso 6")
    newct <- frbind(newct, createEQBCT(newct, bct, stillt))
    rmct <- extractRMCT(constraints(x))
    if (nrow(rmct) != 0L)
      newct <- frbind(newct, rmct)
    x[["constraints"]] <- newct
  }
  x[["parameters"]] <- list(transitions = NULL, states = NULL)
  return(x)
}
