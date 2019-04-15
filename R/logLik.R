# This function returns the logLikelihood

# P : Transition matrix P(will go to i | is in j)
# E : Event observation matrix P(detection event j | is in i)
# e: Sequence of observation events. Since we are taking
# time increase small, it is expected to have mostly missing
# values. The first value is expected to have an observation.

logLik <- function(P, E, e) {
  # Set output vector to zero
  output <- 0
  # Initial steady state calculation
  edecomp <- eigen(P)
  # The steady state might be not unique! Fix this later
  i <- which(abs(edecomp$values - 1) < 1e-6)[1]
  state <- edecomp$vectors[,i]
  state <- state / sum(state)

  # Now run the filter. It is a loop, it should be written in
  # some compiled language for performance.
  for (i in 1:length(e)) {
    if (!is.na(e[i])) {
      state <- E[,e[i]] * state
      output <- output + log(sum(state))
      state <- state / sum(state)
    }
    state <- P %*% matrix(state, ncol = 1)
  }
  return(-output)
}
