#' Adds a transition to the model
#'
#' The specified transition is added to the model as a
#' transition with non zero probability.
#'
#' Since the transition probabilities from the initial state of the
#' newly specified transition still have to sum up to one, the
#' correspondent constraint is modified accordingly.
#'
#' It is not recommended to use this function to define a big model, as
#' it is much slower than specifying all transitions in advance.
#'
#' @param x A HMM object.
#' @param t The transition, as a two dimensional integer vector.
#' The first element is the number of the initial state and the
#' second one the number of the final state.
#'
#' @return A HMM object similar to the input but with the additional
#' transition.
#'
#' @seealso \link{HMM}, \link{addconstraint}
#'
#' @examples
#' model <- HMM(3)
#' model <- addtransition(model, c(1,2))
#' model <- addtransition(model, c(2,3))
#' model <- addtransition(model, c(3,1))
#' transitions(model)
#' constraints(model)
#'
addtransition <- function(x,t) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  if (length(t) < 2)
    stop("Initial and final state must be specified.")
  if (length(t) > 2) {
    warning(paste0("Only initial and final state are required. ",
                   "Remaining elements of the vector are ignored."))
    t <- t[1:2]
  }
  t <- as.integer(t)
  TL <- transitions(x)
  CT <- constraints(x)
  i <- findTorder(TL, t)
  if (i == 0L)
    stop("Transition already stated as having non zero probability.")
  stilli <- which(TL[2, findTorder(TL,c(t[1], 0L)):(findTorder(TL,c(t[1] + 1L, 0L)) - 1)] == t[1])
  stilli <- stilli + findTorder(TL,c(t[1], 0L)) - 2L

  CT <- updateCTaddtran(CT, i - 1L, stilli, length(CT@x),
                        findTorder(TL,c(t[1], 0L)):(findTorder(TL,c(t[1] + 1L, 0L)) - 1L) - 1L)
  newTL <- matrix(0L,nrow = 2, ncol = ncol(TL) + 1)
  newTL[,1:(i - 1)] <- TL[,1:(i - 1)]
  newTL[,i] <- t
  newTL[,(i + 1):(ncol(newTL))] <- TL[,i:(ncol(TL))]
  TL <- newTL

  x[["transitions"]] <- TL
  x[["constraints"]] <- CT
  x[["parameters"]] <- list(transitions = NULL, states = NULL)
  return(x)
}
