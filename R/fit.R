#' Fits a HMM model
#'
#' Fits the transition probabilities of the model by maximum likelihood.
#'
#' The transition probabilities are fitted by ML, subject to the linear constraints
#' specified in the model. It is possible to specify additional non linear constraints,
#' passing the suitable arguments to the optimizer.
#'
#' @param x A HMM model.
#' @param e A vector with the observed events. It admits missing values.
#' @param init Logical specifying whether the initial state found in x is going
#' to be used. Defaults to FALSE, which means that steady state inizialization will be used
#' instead.
#' @param method The optimization algorithm to be used.
#' Defaults to constrOptim from package stats. The other possible choices are donlp2 and
#' solnp.
#' @param ... Arguments to be passed to the optimizer.
#'
#' @return The fitted model.
#'
#' @seealso \link{logLik}, \link{initparams}, \link{minparams}
#'
#' @examples
#' model <- HMMrectangle(20,20)
#' S <- function(x) if (x > 5) return(0) else return(20*log(5/x))
#' emissions(model) <- createEM(c(20,20), towers, S)
#' model <- initparams(model)
#' model <- minparams(model)
#' logLik(model,events)
#' model <- fit(model,events)
#' logLik(model,events)
#'
#' @export
fit <- function(x, e, init = FALSE, method = "constrOptim", ...) {
  if (is.null(x$parameters$reducedparams)) {
    if (is.null(x$parameters$transitions))
      x <- initparams(x)
    x <- minparams(x)
  }
  trmatrix <- createlconmatrix(x$constraints)
  ui <- trmatrix[, -ncol(trmatrix), drop = FALSE]
  ci <- -trmatrix[, ncol(trmatrix)]
  ui <- rbind(ui, -trmatrix[, -ncol(trmatrix), drop = FALSE])
  ci <- c(ci, trmatrix[, ncol(trmatrix)] - rep(1,length(ci)))
  TM <- createTM(x$transitions, x$parameters$transitions, nstates(x))
  if (!init) {
    createsteady(TM)
    ofun <- function(p) return(
      floglik(TM, TM@x, p, x$parameters$reducedparams$transmatrix,
              createsteady(TM, TRUE), emissions(x), as.integer(e) - 1L)
    )
  }
  else
    ofun <- function(p) return(
      floglik(TM, TM@x, p, x$parameters$reducedparams$transmatrix,
              istates(x), emissions(x), as.integer(e) - 1L)
    )
  if (method == "constrOptim") {
    if (!require("stats"))
      stop("Package states required for constrOptim method")
    rparams(x) <- stats::constrOptim(rparams(x), ofun, NULL,
                    ui, ci)$par
  }
  else if (method == "donlp2") {
    if (!require("Rdonlp2"))
      stop("Package Rdonlp2 required for donlp2 method")
    rparams(x) <-
      Rdonlp2::donlp2NLP(rparams(x), ofun,
        par.lower = rep(0, length(rparams(x))),
        par.upper = rep(1, length(rparams(x))),
        ineqA = trmatrix[, -ncol(trmatrix)],
        ineqA.lower = -trmatrix[, ncol(trmatrix)],
        ineqA.upper = 1 - trmatrix[, ncol(trmatrix)], ...)$solution
  }
  else if (method == "solnp") {
    if (!require("Rsolnp"))
      stop("Package Rsolnp required for solnp method")
    rparams(x) <- Rsolnp::solnp(rparams(x), ofun,
                  ineqfun = function(MAT) trmatrix %*% c(MAT, 1),
                  ineqLB = rep(0, nrow(trmatrix)),
                  ineqUB = rep(1, nrow(trmatrix)), ...)$pars
  }
  else stop(paste0("Method unknown: ", method))
if (!init)
  x <- initsteady(x)
return(x)
}
