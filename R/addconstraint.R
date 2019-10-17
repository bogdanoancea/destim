#' Adds a constraint to the model
#'
addconstraint <- function(x, ct) {
  if (class(x) != "HMM")
    stop("This function only works with HMM objects.")
  if (is.matrix(ct)) {
    x[["constraints"]] <- rbind(x[["constraints"]], ct)
  }
  else {
    newct <- numeric(ntransitions(x) + 1)
    newct[ct[1]] <- 1
    newct[ct[2]] <- -1
    x[["constraints"]] <- rbind(x[["constraints"]],
                                matrix(newct, nrow = 1))
  }
  x[["parameters"]] <- list(transitions = NULL, states = NULL)
  return(x)
}
