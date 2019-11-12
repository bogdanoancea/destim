#' A simple Initializer for HMM objects
#'
#' Set initial parameters for an HMM object, as specified.
#'
#' The HMM object contains three fields. States contains information
#' about the states. Transitions contains a list of matrices with
#' two rows. Each column represents a transition with non-zero
#' probability. Transitions in the same matrix have same probability.
#' Emissions is a matrix with the emission probabilities. These are
#' considered fixed.
#'
#' @param x A HMM object.
#'
#' @return A initialized HMM object.
#'
#' @examples
#' HMM(1L, matrix(c(1L,1L), nrow = 2), EM = matrix(1, nrow = 1))
#'
initparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")

  # Sets random parameters before enforcing constraints
  states <- runif(nstates(x))
  transitions <- runif(ntransitions(x))

  # Enforce constraints: we calculate the closest point in the
  # subspace using Lagrange multipliers
  trmatrix <- matrix(as.logical(diag(ntransitions(x))),
                     ncol = ntransitions(x))
  CT <- constraints(x)
  eqcon <- apply(CT, 1, function(irow)
    (irow[ntransitions(x) + 1] == 0) &&
      (sum(irow != 0) == 2) &&
      (sum(irow) == 0))
  EQCT <- CT[eqcon, 1:ntransitions(x), drop = FALSE]
  trmatrixcols <- ncol(trmatrix)
  for (i in seq_len(nrow(EQCT))) {
    variables <- match(c(1, -1), EQCT[i, ], incomparables = 0)
    columns <- c(match(TRUE, trmatrix[variables[1], ]),
                 match(TRUE, trmatrix[variables[2], ]))
    columns <- sort(columns)
    trmatrix[,columns[1]] <- trmatrix[,columns[1]] |
      trmatrix[,columns[2]]
    trmatrix[, columns[2]] <-
      trmatrix[, trmatrixcols, drop = FALSE]
    trmatrixcols <- trmatrixcols - 1
  }
  trmatrix <- trmatrix[,1:trmatrixcols]
  cconstraints <-
    apply(trmatrix, 2, function(v) which(v)[1])
  trmatrix <- matrix(as.numeric(trmatrix), ncol = ncol(trmatrix))
  MCT <- CT[!eqcon, , drop = FALSE]
  MCT <- cbind(MCT[, -ncol(MCT), drop = FALSE] %*% trmatrix,
               MCT[,ncol(MCT), drop = FALSE])
  MCT <- funique(MCT)
  eqsys <- cbind(matrix(0,nrow = nrow(MCT),
                        ncol = nrow(MCT)),
                 MCT)
  eqsys <- rbind(eqsys,
                 cbind(t(MCT[, -ncol(MCT)]),
                     diag(ncol(trmatrix)) * 2,
                     matrix(2 * transitions[cconstraints], ncol = 1)))
  while (TRUE) {
  transitions <-
    solve(eqsys[,1:nrow(eqsys)], eqsys[,ncol(eqsys)])
  transitions <- transitions[-1:(-nrow(MCT))]
  if (all(transitions >= 0 & transitions <= 1))
    break
  transitions[transitions < 0] <- runif(sum(transitions < 0))
  transitions[transitions > 1] <- runif(sum(transitions > 1))
  eqsys[-1:(-nrow(MCT)), ncol(eqsys)] <- 2 * transitions
  }
  states <- states / sum(states)
  x$parameters <- list(states = states,
                       transitions = trmatrix %*% transitions)
  return(x)
}
