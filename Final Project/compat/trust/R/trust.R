trust <- function(objfun, parinit, rinit, rmax, parscale, iterlim = 100,
                  fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps),
                  minimize = TRUE, blather = FALSE, ...) {
  if (missing(parscale)) parscale <- rep(1, length(parinit))
  sgn <- if (minimize) 1 else -1
  ev <- function(x) {
    o <- objfun(x, ...)
    list(value = sgn * o$value, gradient = sgn * as.vector(o$gradient),
         hessian = sgn * o$hessian)
  }
  x <- parinit; r <- rinit
  out <- ev(x)
  if (!is.finite(out$value)) stop("parinit not feasible")
  conv <- FALSE; it <- 0
  for (it in seq_len(iterlim)) {
    g <- out$gradient * parscale
    B <- out$hessian * outer(parscale, parscale)
    e <- eigen((B + t(B)) / 2, symmetric = TRUE)
    lam <- e$values; Q <- e$vectors; gq <- drop(crossprod(Q, g))
    pstep <- function(l) -drop(Q %*% (gq / (lam + l)))
    if (min(lam) > 0 && sqrt(sum(pstep(0)^2)) <= r) {
      p <- pstep(0)
    } else {
      lo <- max(0, -min(lam)) + 1e-10
      h <- function(l) sqrt(sum(pstep(l)^2)) - r
      hi <- lo + 1
      while (h(hi) > 0) hi <- hi * 2 + 1
      if (h(lo) < 0) { p <- pstep(lo); p <- p * r / sqrt(sum(p^2)) }
      else p <- pstep(uniroot(h, c(lo, hi), tol = 1e-12)$root)
    }
    pred <- sum(g * p) + 0.5 * drop(t(p) %*% B %*% p)  # predicted change (negative)
    xnew <- x + p * parscale
    onew <- tryCatch(ev(xnew), error = function(e) list(value = Inf))
    act <- if (is.finite(onew$value)) onew$value - out$value else Inf
    rho <- if (pred < 0) act / pred else -Inf
    if (blather) message(sprintf("iter %d r=%g pred=%g act=%g", it, r, pred, act))
    pnorm2 <- sqrt(sum(p^2))
    if (rho < 0.25) r <- r / 4 else if (rho > 0.75 && pnorm2 > 0.99 * r) r <- min(2 * r, rmax)
    if (is.finite(act) && rho >= 0.1) {
      x <- xnew; out <- onew
      if (abs(pred) < fterm) { conv <- TRUE; break }
    }
    if (abs(pred) < mterm) { conv <- TRUE; break }
  }
  list(value = sgn * out$value, gradient = sgn * out$gradient, hessian = sgn * out$hessian,
       argument = x, converged = conv, iterations = it)
}
