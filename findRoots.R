library(polynom)

polyRoots <- function(y) {
  p <- polynomial(rev(c(1, y)))
  roots <- solve(p)
  return(roots)
}

findRoots <- function(x) {
  alpha <- x[1]
  beta <- x[2]
  gamma <- x[3]
  d <- sqrt(as.complex((alpha ^ 2) - (4 * beta)))
  roots <- c(-gamma, (-alpha + d) / 2, (-alpha - d) / 2)
  return(roots)
}
