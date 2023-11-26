jacob <- function(x) {
  j <- matrix(0, 3, 3)
  j[1, ] <- -1
  j[2, ] <- c(x[2] + x[3], x[1] + x[3], x[1] + x[2])
  j[3, ] <- -c(x[2] * x[3], x[1] * x[3], x[1] * x[2])
  return(j)
}

ffunc <- function(x, y) {
  f <- rep(0, 3)
  f[1] <- x[1] + x[2]  + x[3] + y[1]
  f[2] <- x[1] * x[2] + x[1] * x[3] + x[2] * x[3] - y[2]
  f[3] <- x[1] * x[2] * x[3] + y[3]
  return(f)
}