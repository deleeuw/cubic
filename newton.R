# newton method to solve the equations
# for alpha, beta, gamma

jacobian <- function(x) {
  alpha <- x[1]
  beta <- x[2]
  gamma <- x[3]
  jc <- matrix(c(1, 0, 1, gamma, 1, alpha, 0, gamma, beta), 3, 3, byrow = TRUE)
  return(jc)
}

newton <- function(y,
                   gamma,
                   itmax = 100,
                   eps = 1e-10,
                   verbose = TRUE) {
  p <- y[1]
  q <- y[2]
  r <- y[3]
  itel <- 1
  xold <- c(abFromC(gamma, y), gamma)
  alpha <- xold[1]
  beta <- xold[2]
  gamma <- xold[3]
  f1 <- alpha + gamma - p
  f2 <- beta + (alpha * gamma) - q
  f3 <- (beta * gamma) - r
  fold <- c(f1, f2, f3)
  sold <- lsFunc(xold, y)
  cold <- Inf
  repeat {
    jacginv <- ginv(jacobian(xold))
    delta <- jacginv %*% fold
    xnew <- xold - delta
    cnew <- sqrt(sum(delta ^ 2))
    rate <- cnew / cold
    alpha <- xnew[1]
    beta <- xnew[2]
    gamma <- xnew[3]
    f1 <- alpha + gamma - p
    f2 <- beta + (alpha * gamma) - q
    f3 <- (beta * gamma) - r
    fnew <- c(f1, f2, f3)
    snew <- lsFunc(xnew, y)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, width = 4, format = "d"),
        "fold ",
        formatC(sold, digits = 15, format = "f"),
        "fnew ",
        formatC(snew, digits = 15, format = "f"),
        "cnew ",
        formatC(cnew, digits = 15, format = "f"),
        "rate ",
        formatC(rate, digits = 15, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || (cnew < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    fold <- fnew
    cold <- cnew
  }
  return(list(x = xnew, itel = itel, f = fnew, sigma = snew, rate = rate))
}

newtonRunner <- function(y) {
  gamma0 <- (1:61) / 10 - 3.1
  alpha <- seq(-3, 3, length = 61)
  beta <- seq(-3, 3, length = 61)
  gamma <- seq(-3, 3, length = 61)
  root1 <- seq(-3, 3, length = 61)
  root2 <- seq(-3, 3, length = 61)
  root3 <- seq(-3, 3, length = 61)
  itel <- seq(-3, 3, length = 61)
  sigma <- seq(-3, 3, length = 61)
  rate <- seq(-3, 3, length = 61)
  for (i in 1:61) {
    xtmp <-
      newton(
        y,
        gamma = gamma0[i],
        itmax = 10000,
        verbose = FALSE,
        eps = 1e-10
      )
    roots <- findRoots(xtmp$x)
    alpha[i] <- (xtmp$x)[1]
    beta[i] <- (xtmp$x)[2]
    gamma[i] <- (xtmp$x)[3]
    itel[i] <- xtmp$itel
    sigma[i] <- xtmp$sigma
    rate[i] <- xtmp$rate
    root1[i] <- roots[1]
    root2[i] <- roots[2]
    root3[i] <- roots[3]
  }
  dtmp <- data.frame(
    itel = itel,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    root1 = root1,
    root2 = root2,
    root3 = root3,
    sigma = sigma,
    rate = rate,
    row.names = gamma0
  )
  print(dtmp, digits = 4)
}
