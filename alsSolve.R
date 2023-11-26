# ALS method to solve the equations
# for alpha, beta, gamma

abFromC <- function(gamma, y) {
  p <- y[1]
  q <- y[2]
  r <- y[3]
  u <- matrix(c(1, gamma, 0, 0, 1, gamma), 3, 2)
  v <- c(p - gamma, q, r)
  return(qr.solve(u, v))
}

cFromAb <- function(alpha, beta, y) {
  p <- y[1]
  q <- y[2]
  r <- y[3]
  s <- 1 + (alpha ^ 2) + (beta ^ 2)
  t <- (p - alpha) + alpha * (q - beta) + beta * r
  return(t / s)
}


alsSolve <- function(y,
                     gamma,
                     itmax = 100000,
                     eps = 1e-15,
                     verbose = FALSE) {
  p <- y[1]
  q <- y[2]
  r <- y[3]
  xold <- c(abFromC(gamma, y), gamma)
  fold <- lsFunc(xold, y)
  itel <- 1
  cold <- Inf
  repeat {
    alpha <- xold[1]
    beta <- xold[2]
    gamma <- cFromAb(alpha, beta, y)
    xnew <- c(abFromC(gamma, y), gamma)
    fnew <- lsFunc(xold, y)
    cnew <- sqrt(sum(xold - xnew) ^ 2)
    rate <- cnew / cold
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, digits = 3, format = "d"),
        "fold ",
        formatC(fold, digits = 10, format = "f"),
        "fnew ",
        formatC(fnew, digits = 10, format = "f"),
        "chng ",
        formatC(cnew, digits = 8, format = "f"),
        "rate ",
        formatC(rate, digits = 8, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || (cnew < eps)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    xold <- xnew
    cold <- cnew
  }
  trate <- convRate(anaDer(xnew, y)$hess)
  return(list(
    x = xnew,
    f = fnew,
    rate = rate,
    trate = trate,
    itel = itel
  ))
}

alsRunner <- function(y) {
  gamma0 <- (1:61) / 10 - 3.1
  alpha <- seq(-3, 3, length = 61)
  beta <- seq(-3, 3, length = 61)
  gamma <- seq(-3, 3, length = 61)
  root1 <- seq(-3, 3, length = 61)
  root2 <- seq(-3, 3, length = 61)
  root3 <- seq(-3, 3, length = 61)
  itel <- seq(-3, 3, length = 61)
  rate <- seq(-3, 3, length = 61)
  trate <- seq(-3, 3, length = 61)
  sigma <- seq(-3, 3, length = 61)
  for (i in 1:61) {
    xtmp <-
      alsSolve(
        y,
        gamma = gamma0[i],
        itmax = 10000,
        verbose = FALSE,
        eps = 1e-10
      )
    roots <- findRoots(xtmp$x)
    itel[i] <- xtmp$itel
    rate[i] <- xtmp$rate
    trate[i] <- xtmp$trate
    alpha[i] <- (xtmp$x)[1]
    beta[i] <- (xtmp$x)[2]
    gamma[i] <- (xtmp$x)[3]
    root1[i] <- roots[1]
    root2[i] <- roots[2]
    root3[i] <- roots[3]
    sigma[i] <- xtmp$f
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
    trate = trate,
    row.names = gamma0
  )
  print(dtmp, digits = 4)
}
