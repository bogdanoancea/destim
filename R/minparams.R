#' Minimizes the number of parameters in a HMM model
#'
#' Creates an RPHMM object that reduces the number of parameters
#' of the model to his minimum.
minparams <- function(x) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  trmatrix <- apply(diag(ntransitions(x)), c(1,2), as.logical)
  CT <- getconstraints(x)
  CT <- CT[CT[,ntransitions(x) + 1] == 0,
           1:ntransitions(x), drop = FALSE]
  extrarows <- integer(0)
  for (i in seq_len(nrow(CT))) {
    nonzero <- CT[i, ] != 0
    if (sum(nonzero) == 2 && sum(CT[i, nonzero]) == 0) {
      columns <- c(which(trmatrix[which(nonzero)[1],]),
                   which(trmatrix[which(nonzero)[2],]))
      trmatrix[,columns[1]] <- trmatrix[,columns[1]] |
                                 trmatrix[,columns[2]]
      trmatrix <- trmatrix[,-columns[2], drop = FALSE]
    }
    else extrarows <- c(extrarows, i)
  }
  cconstraints <-
    apply(trmatrix, 2, function(v) which(v)[1])
  trmatrix <- apply(trmatrix, c(1,2), as.numeric)
  MCT <- cbind(CT[extrarows,],rep(0,length(extrarows)))
  CT <- getconstraints(x)
  CT <- CT[CT[,ntransitions(x) + 1] != 0, , drop = FALSE]
  MCT <- rbind(MCT,CT)
  MCT <- cbind(MCT[, -ncol(MCT), drop = FALSE] %*% trmatrix,
               MCT[,ncol(MCT), drop = FALSE])
  qrdecomp <- qr(MCT[, -ncol(MCT), drop = FALSE])
  MCT <- cbind(qr.R(qrdecomp), -t(qr.Q(qrdecomp)) %*% MCT[, ncol(MCT)])
  MCT <- -solve(MCT[,1:nrow(MCT)],MCT[,-(1:nrow(MCT)), drop = FALSE])
  MCT <- rbind(MCT, cbind(diag(length(qrdecomp$pivot) - nrow(MCT)),
                          rep(0,length(qrdecomp$pivot) - nrow(MCT))))
  trmatrix <- trmatrix[, qrdecomp$pivot] %*% MCT
  if (!is.null(x$parameters$reducedparams))
    x$parameters <- x$parameters[c("states", "transitions")]
  x$parameters <- c(x$parameters,list(reducedparams =
                         list(transmatrix = trmatrix,
                              cconstraints = cconstraints,
                              params = NULL)))
  if (!is.null(x$parameters$transitions)) {
    x$parameters$reducedparams$params <- x$parameters$transitions[
      sapply(1:(ncol(trmatrix) - 1), function(x)
        which(apply(trmatrix ** 2, 1, function(y) sum(y) == y[x]))[1])
      ]
    x$parameters$transitions <- (trmatrix %*%
      matrix(c(x$parameters$reducedparams$params, 1), ncol = 1)) [,]
  }
  return(x)
  }
