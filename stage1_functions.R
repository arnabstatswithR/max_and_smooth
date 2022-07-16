c.phi <- 0.8
b.phi <- -1 * c.phi^(-1) * 2^(c.phi - 1) * log(1 - 0.5^c.phi) * (1 - 0.5^c.phi)
a.phi <- -b.phi * log(-log(1 - 0.5^c.phi))

h <- function(xi){a.phi + b.phi * log(-log(1 - (xi + 0.5)^c.phi))}

g <- function(phi){(1 - exp(-exp((phi - a.phi) / b.phi)))^(1/c.phi) - 0.5}

pp.fit.fast <- function(xdat, threshold, npy = 365){
  z <- list()
  n <- length(xdat)
  ex.ind <- xdat > threshold
  n.ex <- sum(ex.ind)
  lrate <- n.ex/n
  xdat.ex <- xdat[ex.ind]
  in2 <- sqrt(6 * var(xdat.ex))/pi
  in1 <- mean(xdat.ex) - 0.57722 * in2
  in3 <- 1e-08
  in2 <- exp(log(in2) + in3 * log(lrate))
  in1 <- threshold - (in2/in3) * (lrate^(-in3) - 1)
  
  log.mu.init <- log(max(in1, exp(-2)))
  log.sigma.by.mu.init <- log(in2) - log.mu.init
  phi.init <- h(0.1)
  
  init <- c(log.mu.init, log.sigma.by.mu.init, phi.init)
  
  pp.lik <- function(a){
    mu <- exp(a[1])
    sc <- exp(a[1] + a[2])
    xi <- g(a[3])
    
    if((1 + xi * (threshold - mu)/sc) < 0){l <- 10^6}else{
      y <- (xdat.ex - mu)/sc
      y <- 1 + xi * y
      if(min(y) <= 0){l <- 10^6}else{
        ll <- log(sc) * n.ex + sum(log(y)) * (1/xi + 1) + 
          n/npy * (1 + xi * (threshold - mu)/sc)^(-1/xi)
        pp <- (4 - c.phi) * log(xi + 0.5) + (4 - 1) * log(0.5 - xi) + 
          (a[3] - a.phi) / b.phi - exp((a[3] - a.phi) / b.phi)
        l <- ll - pp}}
    l}
  
  x <- optim(init, pp.lik, hessian = TRUE, control = list(maxit = 1000))
  
  z$nllh <- x$value
  z$conv <- x$convergence
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z}

#----------------------------------

# MLE

ppfit.mle <- function(station.no, thresholds, Y){
  ppfit.stwise <- pp.fit.fast(Y[station.no, ], 
                              threshold = thresholds[station.no], 
                              npy = 365.25)
  ppfit.stwise}

#----------------------------------

# bootstrapping based covariance matrix of the estimators

rPPP <- function(station.no){
  mu.station <- mles.original[station.no, 1]
  sigma.station <- mles.original[station.no, 2]
  xi.station <- mles.original[station.no, 3]
  threshold.station <- thresholds[station.no]
  tau.station <- sigma.station + xi.station * (threshold.station - mu.station)
  lambda.station <- nt / npy * (1 + xi.station * (threshold.station - mu.station)/
                                  sigma.station)^(-1/xi.station)
  
  log.mu.init <- mles[station.no, 1]
  log.sigma.by.mu.init <- mles[station.no, 2]
  phi.init <- mles[station.no, 3]
  
  init <- c(log.mu.init, log.sigma.by.mu.init, phi.init)
  
  est.PPP <- function(rep.no){
    set.seed(rep.no)
    n.exceed <- rpois(1, lambda.station)
    sample.PPP <- threshold.station + tau.station / xi.station * 
      (runif(n.exceed)^{-xi.station} - 1)
    
    pp.lik <- function(a){
      mu <- exp(a[1])
      sc <- exp(a[1] + a[2])
      xi <- g(a[3])
      
      if((1 + xi * (threshold.station - mu)/sc) < 0){l <- 10^6}else{
        y <- (sample.PPP - mu)/sc
        y <- 1 + xi * y
        if(min(y) <= 0){l <- 10^6}else{
          ll <- log(sc) * n.exceed + sum(log(y)) * (1/xi + 1) + 
            nt/npy * (1 + xi * (threshold.station - mu)/sc)^(-1/xi)
          pp <- (4 - c.phi) * log(xi + 0.5) + (4 - 1) * log(0.5 - xi) + 
            (a[3] - a.phi) / b.phi - exp((a[3] - a.phi) / b.phi)
          l <- ll - pp}}
      l}
    
    x <- optim(init, pp.lik, control = list(maxit = 1000))
    
    if(x$convergence == 0){return(x$par)}else{return(rep(NA, 3))}
  }
  
  ests <- t(sapply(1:1000, est.PPP)) # parametric bootstrap samples = 1000
  final.out <- cov(ests)
  
  final.out}
