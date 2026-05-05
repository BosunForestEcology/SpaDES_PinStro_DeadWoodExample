# Draw next-year DC for each element of DC_vec using a Markov transition matrix.
# transMatrix rows are current DC, columns are next-year DC.
# Row probabilities need not sum to 1; rmultinom normalises.
applyTransition <- function(DC_vec, transMatrix) {
  vapply(DC_vec, function(dc) {
    probs <- transMatrix[dc, ]
    which(stats::rmultinom(1L, 1L, probs) == 1L)
  }, integer(1L))
}
