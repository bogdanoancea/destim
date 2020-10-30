#' @title Number of states.
#'
#' @description Returns the number of states from a HMM object.
#'
#' @param x the HMM object.
#'
#' @return An integer with the number of states of the model.
#'
#' @seealso \link{ntransitions}, \link{nconstraints}
#'
#' @examples
#' model <- HMM(5)
#' nstates(model)
#'
#' @export
nstates <- function(x) {
  UseMethod("nstates")
}
#' @rdname nstates
#' @export
nstates.HMM <- function(x) return(length(x$states$names))

#' @title Number of transitions.
#'
#' @description Returns the number of possible transitions from a HMM object.
#'
#' @param x the HMM object.
#'
#' @return An integer with the number of possible transitions of the model.
#'
#' @seealso \link{transitions}, \link{nstates}, \link{nconstraints}
#'
#' @examples
#' model <- HMM(5)
#' ntransitions(model)
#' model <- addtransition(model, c(1,2))
#' ntransitions(model)
#'
#' @export
ntransitions <- function(x) {
  UseMethod("ntransitions")
}
#' @rdname ntransitions
#' @export
ntransitions.HMM <- function(x) return(ncol(x$transitions))

#' @title Number of constraints.
#'
#' @description Returns the number of constraints from a HMM object.
#'
#' @param x the HMM object.
#'
#' @return An integer with the number of constraints of the model.
#'
#' @seealso \link{constraints}, \link{nstates}, \link{ntransitions}
#'
#' @examples
#' model <- HMM(5)
#' nconstraints(model)
#' model <- HMMrectangle(3, 3)
#' nconstraints(model)
#'
#' @export
nconstraints <- function(x) {
  UseMethod("nconstraints")
}
#' @rdname nconstraints
#' @export
nconstraints.HMM <- function(x) return(nrow(x$constraints))

#' @title Matrix of constraints.
#'
#' @description Returns the matrix of constraints from a HMM object.
#'
#' @param x the HMM object.
#'
#' @return A row major sparse matrix as in \code{\link{HMM}}.
#'
#' @seealso \link{HMM}, \link{nconstraints}, \link{transitions}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' constraints(model)
#' nconstraints(model)
#' nrow(constraints(model)) # should agree
#'
#' @export
constraints <- function(x) {
  UseMethod("constraints")
}
#' @rdname constraints
#' @export
constraints.HMM <- function(x) return(x$constraints)

#' @title Probabilities of transition.
#'
#' @description Returns the probabilities of transition from a HMM object.
#' The object has to be initialized with \code{\link{initparams}},
#' otherwise it will return numeric(0). The order is row major.
#'
#' @param x the HMM object.
#'
#' @return A numeric vector with the probabilities of transition.
#'
#' @seealso \link{HMM}, \link{initparams}, \link{transitions}, \link{ntransitions}
#'
#' @examples
#' model <- HMM(2)
#' model <- addtransition(model,c(1,2))
#' model <- addtransition(model,c(2,1))
#' model <- initparams(model)
#' ptransition(model)
#'
#' @export
ptransition <- function(x) {
  UseMethod("ptransition")
}
#' @rdname ptransition
#' @export
ptransition.HMM <- function(x)
  return(as.numeric(x$parameters$transitions))

#' @title Initial state probabilities.
#'
#' @description Returns the initial state probabilities from a HMM object.
#' The object has to be initialized with \code{\link{initparams}},
#' which generates a random initial state. The vector of probabilities
#' follows the same order as the states, so \code{ptransition(model)[i]} is the
#' probability of state \code{i}. Of course, the probabilities sum up to one.
#'
#' @param x the HMM object.
#'
#' @return A numeric vector with the probabilities.
#'
#' @seealso \link{initparams}, \link{fit}, \link{nstates}, \link{initsteady}
#'
#' @examples
#' model <- HMM(2)
#' model <- addtransition(model,c(1,2))
#' model <- addtransition(model,c(2,1))
#' model <- initparams(model)
#' istates(model)
#' sum(istates(model)) # should be one
#'
#' @export
istates <- function(x) {
  UseMethod("istates")
}
#' @rdname istates
#' @export
istates.HMM <- function(x) return(x$parameters$states)

#' @title Transitions.
#'
#' @description Returns the list of possible transitions from a HMM object.
#' Each column represents a transition, the first row is the initial state and the
#' second row is the final state. The transitions are ordered, first on the initial
#' state and then on the final state. Any transition not listed in the matrix
#' is supposed to be not possible (zero probability).
#'
#' @param x the HMM object.
#'
#' @return An integer matrix with two rows as in \code{\link{HMM}}.
#'
#' @seealso \link{HMM}, \link{ntransitions}, \link{constraints}, \link{ptransition}
#'
#' @examples
#' model <- HMM(2)
#' transitions(model)
#' model <- addtransition(model,c(1,2))
#' model <- addtransition(model,c(2,1))
#' transitions(model)
#' @export
transitions <- function(x) {
  UseMethod("transitions")
}
#' @rdname transitions
#' @export
transitions.HMM <- function(x) return(x$transitions)

#' @title Reduced parameters.
#'
#' @description Returns the values of the minimal set of parameters.
#' The minimal set of parameters are selected by \code{\link{minparams}}.
#' They are a few probabilities of transition that determine the remaining ones
#' because of the constraints. They are used to fit the model.
#'
#' @param x the HMM object.
#'
#' @return A numeric vector with the values of the parameters.
#'
#' @seealso \link{minparams}, \link{ptransition}, \link{gettransmatrix}, \link{fit}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' model <- initparams(model)
#' model <- minparams(model)
#' rparams(model)
#' ntransitions(model)
#' length(rparams(model)) # A much smaller parameter space!
#'
#' @export
rparams <- function(x) {
  UseMethod("rparams")
}
#' @rdname rparams
#' @export
rparams.HMM <- function(x) return(x$parameters$reducedparams$params)

#' @title Transformation matrix
#'
#' @description Returns the transformation matrix that transforms the minimal
#' parameters into the probabilities of transition.
#' The transformation matrix allows obtains the probabilities of transition from
#' the minimal set of parameters. If we append an one at the end of the vector of
#' parameters, the product of this matrix by such vector is the probabilities of
#' transition vector.
#'
#' @param x the HMM object.
#'
#' @return A matrix.
#'
#' @seealso \link{minparams}, \link{ptransition}, \link{rparams}, \link{fit}
#'
#' @examples
#' model <- HMMrectangle(3,3)
#' model <- initparams(model)
#' model <- minparams(model)
#' # Should be close to zero
#' range(ptransition(model) - gettransmatrix(model) %*% c(rparams(model), 1))
#'
#' @export
gettransmatrix <- function(x) {
  UseMethod("gettransmatrix")
}
#' @rdname gettransmatrix
#' @export
gettransmatrix.HMM <- function(x)
  return(x$parameters$reducedparams$transmatrix)

#' @title Emissions matrix
#'
#' @description Returns the matrix of emissions from a HMM object.
#'
#' @param x the HMM object.
#'
#' @return A column major sparse matrix as in \code{\link{HMM}}.
#'
#' @seealso \link{HMM}, \link{emissions<-}, \link{getTM}
#'
#' @examples
#' model <- HMM(2)
#' emissions(model)<-diag(2)
#' emissions(model)
#'
#' @export
emissions <- function(x) {
  UseMethod("emissions")
}
#' @rdname emissions
#' @export
emissions.HMM <- function(x) return(x$emissions)

#' @title Transition matrix.
#'
#' @description Returns the transition matrix from a HMM object.
#'The transition matrix is represented as row major. This way its transpose matrix
#' which is used to left multiply the state is column major.
#'
#' @param x  the HMM object.
#'
#' @return A row major sparse matrix which is the transition matrix of the model.
#'
#' @seealso \link{HMM}, link{emissions}
#'
#' @examples
#' model <- HMM(2)
#' model <- addtransition(model,c(1,2))
#' model <- initparams(model)
#' getTM(model)
#'
#' @export
getTM <- function(x) {
  UseMethod("getTM")
}
#' @rdname getTM
#' @export
getTM.HMM <- function(x) {
  if (is.null(x$parameters$transitions))
     stop(paste0("[destim::getTM] Parameters of the model are ",
      "required to get the transitions matrix"))
  TM <- createTM(x$transitions, x$parameters$transitions, nstates(x))
  return(TM)
}
