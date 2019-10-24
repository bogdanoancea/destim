funique <- function(x) {
i <- 1
j <- nrow(x)
while (i < j) {
  if (min(apply(x[(i + 1):j, , drop = FALSE], 1, function(xrow)
                                     sum((xrow - x[i, ])**2)))
      < sqrt(.Machine$double.eps)) {
    x[i, ] <- x[j, ]
    j <- j - 1
  }
  else {
    i <- i + 1
  }
}
return(x[1:j, ])
}
