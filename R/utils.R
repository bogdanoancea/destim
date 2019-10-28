# These functions are for internal use and won't be exported in the future

funique <- function(x) {
  return(cppfunique(x, sqrt(.Machine$double.eps)))
}
