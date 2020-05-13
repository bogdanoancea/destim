#' This function returns the smooth states

# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first and last value are expected to have
# observations.
# FS: Filtered/predicted states as returned in fstates
#'
#' @export
sstates <- function(...) {
  UseMethod("sstates")
}
#' @rdname sstates
#' @export
sstates.HMM <- function(x, e) {
  fpass <- forward(x,e)
  bpass <- backward(x,e,fpass$scalefactors)
  return(fpass$alpha * bpass)
}
