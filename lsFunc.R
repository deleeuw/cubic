# the least squares loss function
# with its first and second derivatives

lsFunc <- function(x, y) {
  alpha <- x[1]
  beta <- x[2]
  gamma <- x[3]
  p <- y[1]
  q <- y[2]
  r <- y[3]
  f <- 0.0
  f <- f + (alpha + gamma - p) ^ 2
  f <- f + (beta + (alpha * gamma) - q) ^ 2
  f <- f + ((beta * gamma) - r) ^ 2
  return(f / 2)
}

numDer <- function(x, y) {
  g <- grad(lsFunc, x, y = y)
  h <- hessian(lsFunc, x, y = y)
  return(list(grad = g, hess = h))
}

anaDer <- function(x, y) {
  alpha <- x[1]
  beta <- x[2]
  gamma <- x[3]
  p <- y[1]
  q <- y[2]
  r <- y[3]
  da <- (gamma + alpha - p) + gamma * (beta + (alpha * gamma) - q)
  db <- (beta + (alpha * gamma) - q) + gamma * ((beta * gamma) - r)
  dc <-
    (gamma + alpha - p) + alpha * (beta + (alpha * gamma) - q) + beta * ((beta * gamma) - r)
  g <- c(da, db, dc)
  daa <- 1 + gamma ^ 2
  dab <- gamma
  dac <- 1 + (beta - q) + 2 * alpha * gamma
  dbb <- 1 + gamma ^ 2
  dbc <- alpha - r + 2 * beta * gamma
  dcc <- 1 + alpha ^ 2 + beta ^ 2
  h <- matrix(c(daa, dab, dac, dab, dbb, dbc, dac, dbc, dcc), 3, 3)
  return(list(grad = g, hess = h))
}

convRate <- function(hess) {
  a <- hess[3, 3]
  b <- hess[3,][-3]
  c <- hess[1:2, 1:2]
  return(sum(b * solve(c, b)) / a)
}

