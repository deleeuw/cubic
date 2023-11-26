source("findRoots.R")

abFromC <- function(gamma, y) {
  p <- y[1]
  q <- y[2]
  r <- y[3]
  u <- matrix(c(1, gamma, 0, 0, 1, gamma), 3, 2)
  v <- c(p - gamma, q, r)
  return(qr.solve(u, v))
}

newtonLS <-
  function(y,
           gamma,
           itmax = 100,
           eps = 1e-10,
           verbose = FALSE) {
    itel <- 1
    xold <- c(abFromC(gamma, y), gamma)
    sold <- lsFunc(xold, y)
    cold <- Inf
    repeat {
      h <- anaDer(xold, y)
      grad <- h$grad
      hess <- h$hess
      diff <- drop(ginv(hess) %*% grad)
      xnew <- xold - diff
      snew <- lsFunc(xnew, y)
      cnew <- sqrt(sum(diff ^ 2))
      rate <- cnew / cold
      if (verbose)
      {
        cat(
          "itel ", formatC(itel, width = 5, format = "d"),
          "sold ", formatC(sold, digits = 15, width = 20, format = "f"),
          "snew ", formatC(snew, digits = 15, width = 20, format = "f"),
          "cnew ", formatC(cnew, digits = 15, width = 20, format = "f"),
          "rate ", formatC(rate, digits = 15, width = 20, format = "f"),
          "\n")
      }
      if ((itel == itmax) || (cnew < eps)) {
        break
      }
      xold <- xnew
      sold <- snew
      cold <- cnew
      itel <- itel + 1
    }
    return(list(x = xnew, s = snew, itel = itel, rate = rate))
  }

nlsRunner <- function(y) {
  gamma0 <- (1:61) / 10 - 3.1
  alpha <- seq(-3, 3, length = 61)
  beta <- seq(-3, 3, length = 61)
  gamma <- seq(-3, 3, length = 61)
  root1 <- seq(-3, 3, length = 61)
  root2 <- seq(-3, 3, length = 61)
  root3 <- seq(-3, 3, length = 61)
  itel <- seq(-3, 3, length = 61)
  rate <- seq(-3, 3, length = 61)
  sigma <- seq(-3, 3, length = 61)
  for (i in 1:61) {
    xtmp <-
      newtonLS(
        y,
        gamma = gamma0[i],
        itmax = 10000,
        verbose = FALSE,
        eps = 1e-10
      )
    roots <- findRoots(xtmp$x)
    itel[i] <- xtmp$itel
    alpha[i] <- (xtmp$x)[1]
    beta[i] <- (xtmp$x)[2]
    gamma[i] <- (xtmp$x)[3]
    root1[i] <- roots[1]
    root2[i] <- roots[2]
    root3[i] <- roots[3]
    rate[i] <- xtmp$rate
    sigma[i] <- xtmp$s
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
