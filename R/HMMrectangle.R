#' Creates a basic rectangle HMM model
#'
HMMrectangle <- function(x,y) {
  x <- as.integer(x)
  y <- as.integer(y)
  TL <- matrix(0, nrow = 2, ncol = 9*x*y - 6*x - 6*y + 4)
  k <- 1
  for (i in 1:x)
    for (j in 1:y)
      for (newi in (max(i - 1, 1)):(min(i + 1, x)))
        for (newj in (max(j - 1, 1)):(min(j + 1, y))) {
            TL[, k] <- c((j - 1) * x + i, (newj - 1) * x + newi)
            k <- k + 1
        }
  newmodel <- HMM(x*y, TL)
  TL <- TL[2,] - TL[1,]
  CT <- matrix(0, nrow = sum(TL != 0) - length(unique(TL)) + 1,
               ncol = ntransitions(newmodel) + 1)
  k <- 1
  for (i in unique(TL))
    if (i != 0) {
      j <- which(TL == i)
      if (length(j) > 1)
        for (l in 1:(length(j) - 1)) {
          CT[k, j[l]] <- 1
          CT[k, j[l + 1]] <- -1
          k <- k + 1
        }
    }
    newmodel <- addconstraint(newmodel, CT)
    hvmov <- which(TL %in% c(1, -1, x, -x))[1]
    for (i in setdiff(c(1, -1, x, -x), TL[hvmov]))
      if (any(TL == i))
        newmodel <- addconstraint(newmodel, c(hvmov, which(TL == i)[1]))
    diagmov <- which(!(TL %in% c(1,-1, 0, x, -x)))[1]
    if (!is.na(diagmov))
      for (i in setdiff(unique(TL), c(1,-1, 0, x, -x, TL[diagmov])))
        newmodel <- addconstraint(newmodel,
                                  c(diagmov, which(TL == i)[1]))
  return(newmodel)
}
