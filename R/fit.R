#' Fits a HMM model
#'
#' Fits the model by maximum likelihood
#'
fit <- function(x, e, init = FALSE) {
  if (is.null(x$parameters$reducedparams)) {
    if (is.null(x$parameters$transitions))
      x <- initparams(x)
    x <- minparams(x)
  }
  trmatrix <-
    gettransmatrix(x)[x$parameters$reducedparams$cconstraints, ]
  trmatrix <-
    trmatrix[apply(trmatrix,1,function(trow)
      sum(trow[-length(trow)]^2) != 0), ]
  trmatrix <- funique(trmatrix)
  ui <- trmatrix[, -ncol(trmatrix), drop = FALSE]
  ci <- -trmatrix[, ncol(trmatrix)]
  ui <- rbind(ui, -trmatrix[, -ncol(trmatrix), drop = FALSE])
  ci <- c(ci, trmatrix[, ncol(trmatrix)] - rep(1,length(ci)))
  if (!init)
    ofun <- function(p) {
              rparams(x) <- p
              x <- initsteady(x)
              return(logLik(x,e))
            }
  else
    ofun <- function(p) {
      rparams(x) <- p
      return(logLik(x,e))
    }
  rparams(x) <- constrOptim(rparams(x), ofun, NULL, ui, ci)$par
if (!init)
  x <- initsteady(x)
return(x)
}
