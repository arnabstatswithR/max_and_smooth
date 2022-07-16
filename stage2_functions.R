# auxfunctions

mhupdate <- function(acc, att, mh, nattempts = 50, lower = 0.8, higher = 1.2){
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.30) & these.update
  these.high   <- (acc.rate > 0.50) & these.update
  
  mh[these.low]  <- mh[these.low] * lower
  mh[these.high] <- mh[these.high] * higher
  
  acc[these.update] <- 0
  att[these.update] <- 0
  
  results <- list(acc = acc, att = att, mh = mh)
  return(results)
}

#-----------------------------------------------
# update params
#-----------------------------------------------

cormat.update <- function(distmat, rho, nu, fudge = 1e-5){
  cormat <- CorFx(distmat, 1, rho, nu)
  E <- eigen(cormat)
  E$values <- ifelse(E$values < fudge, fudge, E$values)
  cormat <- E$vectors %*% diag(E$values) %*% t(E$vectors)
  chol.cormat <- chol(cormat)
  cormat.logdet <- 2 * sum(log(diag(chol.cormat)))
  cormat.inv <- chol2inv(chol.cormat)
  list(cormat = cormat, chol.cormat = chol.cormat, cormat.inv = cormat.inv, 
       chol.cormat.inv = chol(cormat.inv), cormat.logdet = cormat.logdet)}

cormat.inv.update.inla <- function(rho, c.mat, g1.mat, g2.mat, alpha = 2){
  cormat.inv <- (1 / rho)^4 * c.mat + 2 * (1 / rho)^2 * g1.mat + g2.mat
  tau <- rho^2 / (4 * pi)
  cormat.inv <- tau * cormat.inv
  cormat.logdet <- -2 * sum(log(diag(spam::chol(cormat.inv))))
  list(cormat.inv = cormat.inv, cormat.logdet = cormat.logdet)}

# ns = number of grid cells

# eta.hat ~ Normal(eta, Sigma_eta.hat)
# Q_eta.hat = inverse(Sigma_eta.hat)

# sigmaSq_psi = 1 / exp(logprec_psi)
# sigmaSq_tau = 1 / exp(logprec_tau)
# sigmaSq_phi = 1 / exp(logprec_phi)
# prec_psi = exp(logprec_psi) = 1 / sigmaSq_psi
# prec_tau = exp(logprec_tau) = 1 / sigmaSq_tau
# prec_phi = exp(logprec_phi) = 1 / sigmaSq_phi
# Sigma_eta = bdiag(sigmaSq_psi * diag(ns), sigmaSq_tau * diag(ns), sigmaSq_phi * diag(ns))

# eta ~ Normal(Z %*% nu, Sigma_eta)
# Q_eta = inverse(Sigma_eta)
# Q_eta = bdiag(prec_psi * diag(ns), prec_tau * diag(ns), prec_phi * diag(ns))

# Q_eta.hat_times_eta.hat <- as.vector(Q_eta.hat %*% eta.hat)

# Likelihood: eta.hat ~ Normal(eta, inverse(Q_eta.hat))
# Prior: eta ~ Normal(Z %*% nu, inverse(Q_eta))

eta.update <- function(ns, eta.hat, Q_eta.hat, Q_eta.hat_times_eta.hat, Z.nu,
                       logprec_psi, logprec_tau, logprec_phi){
  
  Q_eta <- bdiag(Diagonal(ns, exp(logprec_psi)), 
                 Diagonal(ns, exp(logprec_tau)), 
                 Diagonal(ns, exp(logprec_phi)))
  
  eta.cov.inv <- Q_eta.hat + Q_eta
  eta.mean.part <- Q_eta.hat_times_eta.hat + as.vector(Q_eta %*% Z.nu)
  
  chol.eta.cov.inv <- spam::chol(eta.cov.inv)
  tchol.eta.cov.inv <- t(chol.eta.cov.inv)
  omega <- spam::forwardsolve(tchol.eta.cov.inv, eta.mean.part)
  mm <- spam::backsolve(chol.eta.cov.inv, omega)
  zz <- rnorm(3 * ns)
  vv <- spam::backsolve(chol.eta.cov.inv, zz)
  eta <- mm + vv
  eta}

# prec_psi = 1 / sigmaSq_psi
# prec_tau = 1 / sigmaSq_tau
# prec_phi = 1 / sigmaSq_phi
# Q_eta = bdiag(prec_psi * diag(ns), prec_tau * diag(ns), prec_phi * diag(ns))

# sd_beta_psi = 100
# sd_beta_tau = 100
# sd_beta_phi = 100

# prec_w_psi = exp(logprec_w_psi)
# prec_w_tau = exp(logprec_w_tau)
# s2_psi = 1 / prec_w_psi
# s2_tau = 1 / prec_w_tau
# rho_psi = exp(logrho_psi)
# rho_tau = exp(logrho_tau)
# Q_rho_psi = Q(rho_psi)
# Q_rho_tau = Q(rho_tau)
# Sigma_nu = bdiag(sd_beta_psi^2, s2_psi * inverse(Q_rho_psi), 
#                  sd_beta_tau^2, s2_tau * inverse(Q_rho_tau),
#                  sd_beta_phi^2)
# Q_nu = bdiag(1 / sd_beta_psi^2, prec_w_psi * Q_rho_psi, 
#              1 / sd_beta_tau^2, prec_w_tau * Q_rho_tau,
#              1 / sd_beta_phi^2)

# Likelihood: eta ~ Normal(Z %*% nu, inverse(Q_eta))
# Prior: nu ~ Normal(0, Q_nu)

nu.update <- function(ns, nmesh, eta, Z, logprec_psi, logprec_tau, logprec_phi, 
                      sd_beta_psi, sd_beta_tau, sd_beta_phi,
                      logprec_w_psi, logprec_w_tau, cormat_psi.inv, cormat_tau.inv){
  
  Q_eta <- bdiag(Diagonal(ns, exp(logprec_psi)), 
                 Diagonal(ns, exp(logprec_tau)), 
                 Diagonal(ns, exp(logprec_phi)))
  
  Q_nu <- bdiag(1 / sd_beta_psi^2, exp(logprec_w_psi) * cormat_psi.inv, 
                1 / sd_beta_tau^2, exp(logprec_w_tau) * cormat_tau.inv,
                1 / sd_beta_phi^2)
  
  tZ_Q_eta <- t(Z) %*% Q_eta
  nu.cov.inv <- Matrix(Q_nu + tZ_Q_eta %*% Z)
  nu.mean.part <- as.vector(tZ_Q_eta %*% eta)
  
  chol.nu.cov.inv <- spam::chol(nu.cov.inv)
  tchol.nu.cov.inv <- t(chol.nu.cov.inv)
  omega <- spam::forwardsolve(tchol.nu.cov.inv, nu.mean.part)
  mm <- spam::backsolve(chol.nu.cov.inv, omega)
  zz <- rnorm(3 + 2 * nmesh)
  vv <- spam::backsolve(chol.nu.cov.inv, zz)
  nu <- mm + vv
  nu}

# Z.nu <- as.vector(Z %*% nu)

# psi <- eta[1:ns]
# tau <- eta[(ns + 1):(2 * ns)]
# phi <- eta[(2 * ns + 1):(3 * ns)]
# Z.nu1 <- Z.nu[1:ns]
# Z.nu2 <- Z.nu[(ns + 1):(2 * ns)]
# Z.nu3 <- Z.nu[(2 * ns + 1):(3 * ns)]

# updated parameter: logprec_psi

# prec_psi = exp(logprec_psi) = 1 / sigmaSq_psi

# Likelihood: psi ~ Normal(Z.nu1, sigmaSq_psi * diag(ns))

# Prior: sigma_psi ~ Exp(lambda_sigma_psi)
# f(logprec_psi) = 0.5 * lambda_sigma_psi * 
#                 exp(-lambda_sigma_psi * exp(-0.5 * logprec_psi) - 0.5 * logprec_psi)
# log[f(logprec_psi)] = c -lambda_sigma_psi * exp(-0.5 * logprec_psi) - 0.5 * logprec_psi

logprec_psi.update <- function(ns, psi, Z.nu1, logprec_psi, lambda_sigma_psi, 
                               att.logprec_psi, acc.logprec_psi, mh.logprec_psi){
  
  att.logprec_psi <- att.logprec_psi + 1
  
  prec_psi <- exp(logprec_psi)
  can.logprec_psi <- logprec_psi + mh.logprec_psi * rnorm(1)
  can.prec_psi <- exp(can.logprec_psi)
  
  can.ll <- 0.5 * ns * can.logprec_psi - 0.5 * can.prec_psi * sum((psi - Z.nu1)^2)
  cur.ll <- 0.5 * ns * logprec_psi - 0.5 * prec_psi * sum((psi - Z.nu1)^2)
  
  can.lprior <- -lambda_sigma_psi * exp(-0.5 * can.logprec_psi) - 0.5 * can.logprec_psi
  cur.lprior <- -lambda_sigma_psi * exp(-0.5 * logprec_psi) - 0.5 * logprec_psi
  
  ratio <- can.ll + can.lprior - cur.ll - cur.lprior
  
  if(log(runif(1)) < ratio){
    logprec_psi <- can.logprec_psi
    acc.logprec_psi <- acc.logprec_psi + 1}
  
  results <- list(logprec_psi = logprec_psi,
                  att.logprec_psi = att.logprec_psi, 
                  acc.logprec_psi = acc.logprec_psi)
  results}

# updated parameter: logprec_tau

# prec_tau = exp(logprec_tau) = 1 / sigmaSq_tau

# Likelihood: tau ~ Normal(Z.nu2, sigmaSq_tau * diag(ns))

# Prior: sigma_tau ~ Exp(lambda_sigma_tau)
# f(logprec_tau) = 0.5 * lambda_sigma_tau * 
#                 exp(-lambda_sigma_tau * exp(-0.5 * logprec_tau) - 0.5 * logprec_tau)
# log[f(logprec_tau)] = c -lambda_sigma_tau * exp(-0.5 * logprec_tau) - 0.5 * logprec_tau

logprec_tau.update <- function(ns, tau, Z.nu2, logprec_tau, lambda_sigma_tau, 
                               att.logprec_tau, acc.logprec_tau, mh.logprec_tau){
  
  att.logprec_tau <- att.logprec_tau + 1
  
  prec_tau <- exp(logprec_tau)
  can.logprec_tau <- logprec_tau + mh.logprec_tau * rnorm(1)
  can.prec_tau <- exp(can.logprec_tau)
  
  can.ll <- 0.5 * ns * can.logprec_tau - 0.5 * can.prec_tau * sum((tau - Z.nu2)^2)
  cur.ll <- 0.5 * ns * logprec_tau - 0.5 * prec_tau * sum((tau - Z.nu2)^2)
  
  can.lprior <- -lambda_sigma_tau * exp(-0.5 * can.logprec_tau) - 0.5 * can.logprec_tau
  cur.lprior <- -lambda_sigma_tau * exp(-0.5 * logprec_tau) - 0.5 * logprec_tau
  
  ratio <- can.ll + can.lprior - cur.ll - cur.lprior
  
  if(log(runif(1)) < ratio){
    logprec_tau <- can.logprec_tau
    acc.logprec_tau <- acc.logprec_tau + 1}
  
  results <- list(logprec_tau = logprec_tau,
                  att.logprec_tau = att.logprec_tau, 
                  acc.logprec_tau = acc.logprec_tau)
  results}

# updated parameter: logprec_phi

# prec_phi = exp(logprec_phi) = 1 / sigmaSq_phi

# Likelihood: phi ~ Normal(Z.nu3, sigmaSq_phi * diag(ns))

# Prior: sigma_phi ~ Exp(lambda_sigma_phi)
# f(logprec_phi) = 0.5 * lambda_sigma_phi * 
#                 exp(-lambda_sigma_phi * exp(-0.5 * logprec_phi) - 0.5 * logprec_phi)
# log[f(logprec_phi)] = c -lambda_sigma_phi * exp(-0.5 * logprec_phi) - 0.5 * logprec_phi

logprec_phi.update <- function(ns, phi, Z.nu3, logprec_phi, lambda_sigma_phi, 
                               att.logprec_phi, acc.logprec_phi, mh.logprec_phi){
  
  att.logprec_phi <- att.logprec_phi + 1
  
  prec_phi <- exp(logprec_phi)
  can.logprec_phi <- logprec_phi + mh.logprec_phi * rnorm(1)
  can.prec_phi <- exp(can.logprec_phi)
  
  can.ll <- 0.5 * ns * can.logprec_phi - 0.5 * can.prec_phi * sum((phi - Z.nu3)^2)
  cur.ll <- 0.5 * ns * logprec_phi - 0.5 * prec_phi * sum((phi - Z.nu3)^2)
  
  can.lprior <- -lambda_sigma_phi * exp(-0.5 * can.logprec_phi) - 0.5 * can.logprec_phi
  cur.lprior <- -lambda_sigma_phi * exp(-0.5 * logprec_phi) - 0.5 * logprec_phi
  
  ratio <- can.ll + can.lprior - cur.ll - cur.lprior
  
  if(log(runif(1)) < ratio){
    logprec_phi <- can.logprec_phi
    acc.logprec_phi <- acc.logprec_phi + 1}
  
  results <- list(logprec_phi = logprec_phi,
                  att.logprec_phi = att.logprec_phi, 
                  acc.logprec_phi = acc.logprec_phi)
  results}

# w_psi <- nu[2:(nmesh + 1)]
# w_tau <- nu[(nmesh + 3):(2 * nmesh + 2)]

# ss.w_psi <- sum(w_psi * as.vector(cormat_psi.inv %*% w_psi))
# ss.w_tau <- sum(w_tau * as.vector(cormat_tau.inv %*% w_tau))

# updated parameter: logprec_w_psi, logrho_psi

# prec_w_psi = exp(logprec_w_psi)
# s2_psi = 1 / prec_w_psi
# rho_psi = exp(logrho_psi)

# Likelihood: w_psi ~ Normal(0, s2_psi * inverse(cormat_psi.inv))

# cur.ss.w_psi = sum(w_psi * as.vector(cur.cormat_psi.inv %*% w_psi))

# Prior: s_psi ~ Exp(lambda_s_psi), rho_psi ~ Inverse-gamma(1, lambda_rho_psi)
# f(logprec_w_psi) = 0.5 * lambda_s_psi * 
#                    exp(-lambda_s_psi * exp(-0.5 * logprec_w_psi) - 0.5 * logprec_w_psi)
# log[f(logprec_w_psi)] = c -lambda_s_psi * exp(-0.5 * logprec_w_psi) - 0.5 * logprec_w_psi

# f(logrho_psi) = lambda_rho_psi * 
#                    exp(-lambda_rho_psi * exp(-logrho_psi) - logrho_psi)
# log[f(logrho_psi)] = c -lambda_rho_psi * exp(-logrho_psi) - logrho_psi

# cur.cormat_psi.details <- cormat.inv.update.inla(rho_psi, c.mat, g1.mat, g2.mat, alpha = alpha)
# cur.cormat_psi.inv <- cur.cormat_psi.details$cormat.inv
# cur.cormat_psi.logdet <- cur.cormat_psi.details$cormat.logdet

logprec_w_psi_logrho_psi.update <- function(nmesh, logprec_w_psi, logrho_psi, w_psi, 
                                            lambda_s_psi, lambda_rho_psi,
                                            cur.cormat_psi.inv, cur.cormat_psi.logdet,
                                            cur.ss.w_psi, c.mat, g1.mat, g2.mat,
                                            att.logprec_w_psi_logrho_psi, 
                                            acc.logprec_w_psi_logrho_psi, 
                                            mh.logprec_w_psi_logrho_psi){
  
  att.logprec_w_psi_logrho_psi <- att.logprec_w_psi_logrho_psi + 1
  
  prec_w_psi <- exp(logprec_w_psi)
  can.logprec_w_psi <- logprec_w_psi + mh.logprec_w_psi_logrho_psi * rnorm(1)
  can.prec_w_psi <- exp(can.logprec_w_psi)
  
  rho_psi <- exp(logrho_psi)
  can.logrho_psi <- logrho_psi + mh.logprec_w_psi_logrho_psi * rnorm(1)
  can.rho_psi <- exp(can.logrho_psi)
  
  can.cormat_psi.details <- cormat.inv.update.inla(can.rho_psi, c.mat, g1.mat, g2.mat)
  
  can.cormat_psi.inv <- can.cormat_psi.details$cormat.inv
  can.cormat_psi.logdet <- can.cormat_psi.details$cormat.logdet
  
  can.ss.w_psi <- sum(w_psi * as.vector(can.cormat_psi.inv %*% w_psi))
  
  can.ll <- 0.5 * nmesh * can.logprec_w_psi - 
    0.5 * can.cormat_psi.logdet - 0.5 * can.prec_w_psi * can.ss.w_psi
  cur.ll <- 0.5 * nmesh * logprec_w_psi - 
    0.5 * cur.cormat_psi.logdet - 0.5 * prec_w_psi * cur.ss.w_psi
  
  can.lprior <- -lambda_s_psi * exp(-0.5 * can.logprec_w_psi) - 0.5 * can.logprec_w_psi - 
    lambda_rho_psi * exp(-can.logrho_psi) - can.logrho_psi
  cur.lprior <- -lambda_s_psi * exp(-0.5 * logprec_w_psi) - 0.5 * logprec_w_psi - 
    lambda_rho_psi * exp(-logrho_psi) - logrho_psi
  
  ratio <- can.ll + can.lprior - cur.ll - cur.lprior
  
  if(log(runif(1)) < ratio){
    logprec_w_psi <- can.logprec_w_psi
    logrho_psi <- can.logrho_psi
    cur.cormat_psi.inv <- can.cormat_psi.inv
    cur.cormat_psi.logdet <- can.cormat_psi.logdet
    cur.ss.w_psi <- can.ss.w_psi
    acc.logprec_w_psi_logrho_psi <- acc.logprec_w_psi_logrho_psi + 1}
  
  results <- list(logprec_w_psi = logprec_w_psi,
                  logrho_psi = logrho_psi,
                  cur.cormat_psi.inv = cur.cormat_psi.inv,
                  cur.cormat_psi.logdet = cur.cormat_psi.logdet,
                  cur.ss.w_psi = cur.ss.w_psi,
                  att.logprec_w_psi_logrho_psi = att.logprec_w_psi_logrho_psi, 
                  acc.logprec_w_psi_logrho_psi = acc.logprec_w_psi_logrho_psi)
  results}

# updated parameter: logprec_w_tau, logrho_tau

# prec_w_tau = exp(logprec_w_tau)
# s2_tau = 1 / prec_w_tau
# rho_tau = exp(logrho_tau)

# Likelihood: w_tau ~ Normal(0, s2_tau * inverse(cormat_tau.inv))

# cur.ss.w_tau = sum(w_tau * as.vector(cur.cormat_tau.inv %*% w_tau))

# Prior: s_tau ~ Exp(lambda_s_tau), rho_tau ~ Inverse-gamma(1, lambda_rho_tau)
# f(logprec_w_tau) = 0.5 * lambda_s_tau * 
#                    exp(-lambda_s_tau * exp(-0.5 * logprec_w_tau) - 0.5 * logprec_w_tau)
# log[f(logprec_w_tau)] = c -lambda_s_tau * exp(-0.5 * logprec_w_tau) - 0.5 * logprec_w_tau

# f(logrho_tau) = lambda_rho_tau * 
#                    exp(-lambda_rho_tau * exp(-logrho_tau) - logrho_tau)
# log[f(logrho_tau)] = c -lambda_rho_tau * exp(-logrho_tau) - logrho_tau

# cur.cormat_tau.details <- cormat.inv.update.inla(rho_tau, c.mat, g1.mat, g2.mat, alpha = alpha)
# cur.cormat_tau.inv <- cur.cormat_tau.details$cormat.inv
# cur.cormat_tau.logdet <- cur.cormat_tau.details$cormat.logdet

logprec_w_tau_logrho_tau.update <- function(nmesh, logprec_w_tau, logrho_tau, w_tau, 
                                            lambda_s_tau, lambda_rho_tau,
                                            cur.cormat_tau.inv, cur.cormat_tau.logdet,
                                            cur.ss.w_tau, c.mat, g1.mat, g2.mat,
                                            att.logprec_w_tau_logrho_tau, 
                                            acc.logprec_w_tau_logrho_tau, 
                                            mh.logprec_w_tau_logrho_tau){
  
  att.logprec_w_tau_logrho_tau <- att.logprec_w_tau_logrho_tau + 1
  
  prec_w_tau <- exp(logprec_w_tau)
  can.logprec_w_tau <- logprec_w_tau + mh.logprec_w_tau_logrho_tau * rnorm(1)
  can.prec_w_tau <- exp(can.logprec_w_tau)
  
  rho_tau <- exp(logrho_tau)
  can.logrho_tau <- logrho_tau + mh.logprec_w_tau_logrho_tau * rnorm(1)
  can.rho_tau <- exp(can.logrho_tau)
  
  can.cormat_tau.details <- cormat.inv.update.inla(can.rho_tau, c.mat, g1.mat, g2.mat)
  
  can.cormat_tau.inv <- can.cormat_tau.details$cormat.inv
  can.cormat_tau.logdet <- can.cormat_tau.details$cormat.logdet
  
  can.ss.w_tau <- sum(w_tau * as.vector(can.cormat_tau.inv %*% w_tau))
  
  can.ll <- 0.5 * nmesh * can.logprec_w_tau - 
    0.5 * can.cormat_tau.logdet - 0.5 * can.prec_w_tau * can.ss.w_tau
  cur.ll <- 0.5 * nmesh * logprec_w_tau - 
    0.5 * cur.cormat_tau.logdet - 0.5 * prec_w_tau * cur.ss.w_tau
  
  can.lprior <- -lambda_s_tau * exp(-0.5 * can.logprec_w_tau) - 0.5 * can.logprec_w_tau - 
    lambda_rho_tau * exp(-can.logrho_tau) - can.logrho_tau
  cur.lprior <- -lambda_s_tau * exp(-0.5 * logprec_w_tau) - 0.5 * logprec_w_tau - 
    lambda_rho_tau * exp(-logrho_tau) - logrho_tau
  
  ratio <- can.ll + can.lprior - cur.ll - cur.lprior
  
  if(log(runif(1)) < ratio){
    logprec_w_tau <- can.logprec_w_tau
    logrho_tau <- can.logrho_tau
    cur.cormat_tau.inv <- can.cormat_tau.inv
    cur.cormat_tau.logdet <- can.cormat_tau.logdet
    cur.ss.w_tau <- can.ss.w_tau
    acc.logprec_w_tau_logrho_tau <- acc.logprec_w_tau_logrho_tau + 1}
  
  results <- list(logprec_w_tau = logprec_w_tau,
                  logrho_tau = logrho_tau,
                  cur.cormat_tau.inv = cur.cormat_tau.inv,
                  cur.cormat_tau.logdet = cur.cormat_tau.logdet,
                  cur.ss.w_tau = cur.ss.w_tau,
                  att.logprec_w_tau_logrho_tau = att.logprec_w_tau_logrho_tau, 
                  acc.logprec_w_tau_logrho_tau = acc.logprec_w_tau_logrho_tau)
  results}

#-----------------------------------------------
# mcmc
#-----------------------------------------------

# inputs

# mles.covmats: A list of length same as the number of spatial locations,
#               i-th element includes the MLEs (vector of length 3) of 
#               log-mu, log-sigma/mu, phi from Stage 1, and also includes 
#               the covariance matrix (3-by-3 dimensional) of the estimators

# loc: the ns-times-2 matrix where each row represents a coordinate of a grid cell

# inla.mats: it includes the matrices c.mat, g1.mat, g2.mat, A returned by 
#             inla.mesh.fem and inla.spde.make.A functions from INLA

# alpha: alpha = 2 indicates we set the Matern smoothness to (alpha - 1) = 1

# beta_psi.init: initial value for beta_psi
# beta_tau.init: initial value for beta_tau
# beta_phi.init: initial value of beta_psi
# w_psi.init: initial value of the vector w_psi (denoted by w_psi_star in the chapter)
# w_tau.init: initial value of the vector w_tau (denoted by w_tau_star in the chapter)
# sigmaSq_psi.init: initial value of sigmaSq_psi
# sigmaSq_tau.init: initial value of sigmaSq_tau
# sigmaSq_phi.init: initial value of sigmaSq_phi
# s2_psi.init: initial value of s2_psi (denoted by s^2_phi in the chapter)
# s2_tau.init: initial value of s2_tau (denoted by s^2_tau in the chapter)
# rho_psi.init: initial value of rho_psi
# rho_tau.init: initial value of rho_tau

# prior choices

# sd_beta_psi: prior standard deviation of beta_psi, set to 100
# sd_beta_tau: prior standard deviation of beta_tau, set to 100
# sd_beta_phi: prior standard deviation of beta_phi, set to 100
# lambda_sigma_psi: prior rate parameter for sigma_psi ~ Exp(lambda_sigma_psi)
# lambda_sigma_tau: prior rate parameter for sigma_tau ~ Exp(lambda_sigma_tau)
# lambda_sigma_phi: prior rate parameter for sigma_phi ~ Exp(lambda_sigma_phi)
# lambda_s_psi: prior rate parameter for s_psi ~ Exp(lambda_s_psi)
# lambda_s_tau: prior rate parameter for s_tau ~ Exp(lambda_s_tau)
# lambda_rho_psi: prior rate parameter for rho_psi ~ Inverse-gamma(1, lambda_rho_psi)
# lambda_rho_tau: prior rate parameter for rho_tau ~ Inverse-gamma(1, lambda_rho_tau)

# outputs: MCMC chains from all the parameters, smooth parameter profiles for selected stations,
#          posterior mean and variance of smooth parameter profiles for all stations,
#          posterior mean and variance of w_psi and w_tau

mcmc.maxNsmooth <- function(mles.covmats, loc, inla.mats, alpha = 2,
                            beta_psi.init = NULL, beta_tau.init = NULL,
                            beta_phi.init = NULL, w_psi.init = NULL,
                            w_tau.init = NULL, sigmaSq_psi.init = NULL,
                            sigmaSq_tau.init = NULL, sigmaSq_phi.init = NULL,
                            s2_psi.init = NULL, s2_tau.init = NULL,
                            rho_psi.init = NULL, rho_tau.init = NULL,
                            # priors
                            sd_beta_psi = 1e2, 
                            sd_beta_tau = 1e2, 
                            sd_beta_phi = 1e2,
                            lambda_sigma_psi = 0.1,
                            lambda_sigma_tau = 0.1,
                            lambda_sigma_phi = 0.1,
                            lambda_s_psi = 0.1,
                            lambda_s_tau = 0.1,
                            lambda_rho_psi = 0.1,
                            lambda_rho_tau = 0.1,
                            # mcmc settings
                            iters = 4000, burn = 2000, thin = 5){
  
  tick <- proc.time()[3]
  
  library(Matrix)
  library(spam)
  library(fields)
  
  c.mat <- inla.mats$c.mat
  g1.mat <- inla.mats$g1.mat
  g2.mat <- inla.mats$g2.mat
  A <- inla.mats$A
  
  ns <- nrow(loc)
  nmesh <- ncol(A)
  
  # Make the Z (design) matrix
  
  Z <- bdiag(cbind(1, A), cbind(1, A), rep(1, ns))
  
  # extract necessary elements from mles.covmats
  
  psi <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$mle[1]})
  tau <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$mle[2]})
  phi <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$mle[3]})
  v_psi <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$covmat[1, 1]})
  v_tau <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$covmat[2, 2]})
  v_phi <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$covmat[3, 3]})
  v_psi_tau <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$covmat[1, 2]})
  v_psi_phi <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$covmat[1, 3]})
  v_tau_phi <- sapply(1:ns, function(loc.no){mles.covmats[[loc.no]]$covmat[2, 3]})
  
  # make the vector and covariance matrix for the observations
  
  eta.hat <- c(psi, tau, phi)
  
  library(Matrix)
  
  Sigma_eta.hat <- matrix(0, nrow = 3 * ns, ncol = 3 * ns)
  for(i in 1:ns){
    Sigma_eta.hat[i, i] <- v_psi[i]
    Sigma_eta.hat[i + ns, i + ns] <- v_tau[i]
    Sigma_eta.hat[i + 2 * ns, i + 2 * ns] <- v_phi[i]
    
    Sigma_eta.hat[i, i + ns] <- Sigma_eta.hat[i + ns, i] <- v_psi_tau[i]
    Sigma_eta.hat[i, i + 2 * ns] <- Sigma_eta.hat[i + 2 * ns, i] <- v_psi_phi[i]
    Sigma_eta.hat[i + ns, i + 2 * ns] <- Sigma_eta.hat[i + 2 * ns, i + ns] <- v_tau_phi[i]
  }
  Sigma_eta.hat <- Matrix(Sigma_eta.hat)
  Q_eta.hat <- solve(Sigma_eta.hat)
  
  Q_eta.hat_times_eta.hat <- as.vector(Q_eta.hat %*% eta.hat)
  
  if(is.null(beta_psi.init)){
    beta_psi <- sum(psi / v_psi) / sum(1 / v_psi)
  }else{beta_psi <- beta_psi.init}
  
  if(is.null(beta_tau.init)){
    beta_tau <- sum(tau / v_tau) / sum(1 / v_tau)
  }else{beta_tau <- beta_tau.init}
  
  if(is.null(beta_phi.init)){
    beta_phi <- sum(phi / v_phi) / sum(1 / v_phi)
  }else{beta_phi <- beta_phi.init}
  
  eta <- c(rep(beta_psi, ns), rep(beta_tau, ns), rep(beta_phi, ns))
  
  if(is.null(sigmaSq_psi.init)){
    sigmaSq_psi <- var(psi - beta_psi) / 2
  }else{sigmaSq_psi <- sigmaSq_psi.init}
  
  if(is.null(sigmaSq_tau.init)){
    sigmaSq_tau <- var(tau - beta_tau) / 2
  }else{sigmaSq_tau <- sigmaSq_tau.init}
  
  if(is.null(sigmaSq_phi.init)){
    sigmaSq_phi <- var(phi - beta_phi) / 2
  }else{sigmaSq_phi <- sigmaSq_phi.init}
  
  if(is.null(w_psi.init)){
    w_psi <- rep(0, nmesh)
  }else{w_psi <- w_psi.init}
  
  if(is.null(w_tau.init)){
    w_tau <- rep(0, nmesh)
  }else{w_tau <- w_tau.init}
  
  if(is.null(s2_psi.init)){
    s2_psi <- var(psi - beta_psi) / 2
  }else{s2_psi <- s2_psi.init}
  
  if(is.null(s2_tau.init)){
    s2_tau <- var(tau - beta_tau) / 2
  }else{s2_tau <- s2_tau.init}
  
  nu <- c(beta_psi, w_psi, beta_tau, w_tau, beta_phi)
  
  Z.nu <- as.vector(Z %*% nu)
  
  psi <- eta[1:ns]
  tau <- eta[(ns + 1):(2 * ns)]
  phi <- eta[(2 * ns + 1):(3 * ns)]
  Z.nu1 <- Z.nu[1:ns]
  Z.nu2 <- Z.nu[(ns + 1):(2 * ns)]
  Z.nu3 <- Z.nu[(2 * ns + 1):(3 * ns)]
  
  max.dist <- max(rdist(loc)) # maximum distance between two points
  
  if(is.null(rho_psi.init)){
    rho_psi <- 0.05 * max.dist
  }else{rho_psi <- rho_psi.init}
  
  if(is.null(rho_tau.init)){
    rho_tau <- 0.05 * max.dist
  }else{rho_tau <- rho_tau.init}
  
  cur.cormat_psi.details <- cormat.inv.update.inla(rho_psi, c.mat, g1.mat, g2.mat, alpha = alpha)
  cur.cormat_psi.inv <- cur.cormat_psi.details$cormat.inv
  cur.cormat_psi.logdet <- cur.cormat_psi.details$cormat.logdet
  
  cur.cormat_tau.details <- cormat.inv.update.inla(rho_tau, c.mat, g1.mat, g2.mat, alpha = alpha)
  cur.cormat_tau.inv <- cur.cormat_tau.details$cormat.inv
  cur.cormat_tau.logdet <- cur.cormat_tau.details$cormat.logdet
  
  beta_psi <- nu[1] 
  w_psi <- nu[2:(nmesh + 1)] 
  beta_tau <- nu[nmesh + 2]
  w_tau <- nu[(nmesh + 3):(2 * nmesh + 2)]
  beta_phi <- nu[2 * nmesh + 3]
  
  cur.ss.w_psi <- sum(w_psi * as.vector(cur.cormat_psi.inv %*% w_psi))
  cur.ss.w_tau <- sum(w_tau * as.vector(cur.cormat_tau.inv %*% w_tau))
  
  logprec_psi <- -log(sigmaSq_psi)
  logprec_tau <- -log(sigmaSq_tau)
  logprec_phi <- -log(sigmaSq_phi)
  logprec_w_psi <- -log(s2_psi)
  logprec_w_tau <- -log(s2_tau)
  logrho_psi <- log(rho_psi)
  logrho_tau <- log(rho_tau)
  
  att.logprec_psi <- acc.logprec_psi <- mh.logprec_psi <- 0.1
  att.logprec_tau <- acc.logprec_tau <- mh.logprec_tau <- 0.1
  att.logprec_phi <- acc.logprec_phi <- mh.logprec_phi <- 0.1
  att.logprec_w_psi_logrho_psi <- acc.logprec_w_psi_logrho_psi <- mh.logprec_w_psi_logrho_psi <- 0.1
  att.logprec_w_tau_logrho_tau <- acc.logprec_w_tau_logrho_tau <- mh.logprec_w_tau_logrho_tau <- 0.1
  
  keepers.beta_psi <- rep(NA, iters)
  keepers.beta_tau <- rep(NA, iters)
  keepers.beta_phi <- rep(NA, iters)
  keepers.sigma_psi <- rep(NA, iters)
  keepers.sigma_tau <- rep(NA, iters)
  keepers.sigma_phi <- rep(NA, iters)
  keepers.s_psi <- rep(NA, iters)
  keepers.s_tau <- rep(NA, iters)
  keepers.rho_psi <- rep(NA, iters)
  keepers.rho_tau <- rep(NA, iters)
  
  keepers.psi <- matrix(NA, iters, ns)
  keepers.tau <- matrix(NA, iters, ns)
  keepers.phi <- matrix(NA, iters, ns)
  
  psi.sum <- psi2.sum <- rep(0, ns)
  tau.sum <- tau2.sum <- rep(0, ns)
  phi.sum <- phi2.sum <- rep(0, ns)
  
  w_psi.sum <- w_psi2.sum <- rep(0, nmesh)
  w_tau.sum <- w_tau2.sum <- rep(0, nmesh)
  
  return.iters <- (burn + 1):iters
  
  for(iter in 1:iters){for(ttt in 1:thin){
    
    # updated parameter: eta
    
    eta <- eta.update(ns, eta.hat, Q_eta.hat, Q_eta.hat_times_eta.hat, Z.nu,
                      logprec_psi, logprec_tau, logprec_phi)
    
    psi <- eta[1:ns]
    tau <- eta[(ns + 1):(2 * ns)]
    phi <- eta[(2 * ns + 1):(3 * ns)]
    
    # updated parameter: nu
    
    nu <- nu.update(ns, nmesh, eta, Z, logprec_psi, logprec_tau, logprec_phi, 
                    sd_beta_psi, sd_beta_tau, sd_beta_phi,
                    logprec_w_psi, logprec_w_tau, cur.cormat_psi.inv, cur.cormat_tau.inv)
    
    # some calculations for further updates
    
    Z.nu <- as.vector(Z %*% nu)
    
    Z.nu1 <- Z.nu[1:ns]
    Z.nu2 <- Z.nu[(ns + 1):(2 * ns)]
    Z.nu3 <- Z.nu[(2 * ns + 1):(3 * ns)]
    
    beta_psi <- nu[1] 
    w_psi <- nu[2:(nmesh + 1)] 
    beta_tau <- nu[nmesh + 2]
    w_tau <- nu[(nmesh + 3):(2 * nmesh + 2)]
    beta_phi <- nu[2 * nmesh + 3]
    
    cur.ss.w_psi <- sum(w_psi * as.vector(cur.cormat_psi.inv %*% w_psi))
    cur.ss.w_tau <- sum(w_tau * as.vector(cur.cormat_tau.inv %*% w_tau))
    
    # updated parameter: logprec_psi
    
    logprec_psi.update.details <- 
      logprec_psi.update(ns, psi, Z.nu1, logprec_psi, lambda_sigma_psi, 
                         att.logprec_psi, acc.logprec_psi, mh.logprec_psi)
    
    logprec_psi <- logprec_psi.update.details$logprec_psi
    att.logprec_psi <- logprec_psi.update.details$att.logprec_psi
    acc.logprec_psi <- logprec_psi.update.details$acc.logprec_psi
    
    prec_psi <- exp(logprec_psi)
    sigmaSq_psi <- 1 / prec_psi
    sigma_psi <- sqrt(sigmaSq_psi)
    
    # updated parameter: logprec_tau
    
    logprec_tau.update.details <- 
      logprec_tau.update(ns, tau, Z.nu2, logprec_tau, lambda_sigma_tau, 
                         att.logprec_tau, acc.logprec_tau, mh.logprec_tau)
    
    logprec_tau <- logprec_tau.update.details$logprec_tau
    att.logprec_tau <- logprec_tau.update.details$att.logprec_tau
    acc.logprec_tau <- logprec_tau.update.details$acc.logprec_tau
    
    prec_tau <- exp(logprec_tau)
    sigmaSq_tau <- 1 / prec_tau
    sigma_tau <- sqrt(sigmaSq_tau)
    
    # updated parameter: logprec_phi
    
    logprec_phi.update.details <- 
      logprec_phi.update(ns, phi, Z.nu3, logprec_phi, lambda_sigma_phi, 
                         att.logprec_phi, acc.logprec_phi, mh.logprec_phi)
    
    logprec_phi <- logprec_phi.update.details$logprec_phi
    att.logprec_phi <- logprec_phi.update.details$att.logprec_phi
    acc.logprec_phi <- logprec_phi.update.details$acc.logprec_phi
    
    prec_phi <- exp(logprec_phi)
    sigmaSq_phi <- 1 / prec_phi
    sigma_phi <- sqrt(sigmaSq_phi)
    
    # updated parameter: logprec_w_psi, logrho_psi
    
    logprec_w_psi_logrho_psi.update.details <- 
      logprec_w_psi_logrho_psi.update(nmesh, logprec_w_psi, logrho_psi, w_psi, 
                                      lambda_s_psi, lambda_rho_psi,
                                      cur.cormat_psi.inv, cur.cormat_psi.logdet,
                                      cur.ss.w_psi, c.mat, g1.mat, g2.mat,
                                      att.logprec_w_psi_logrho_psi, 
                                      acc.logprec_w_psi_logrho_psi, 
                                      mh.logprec_w_psi_logrho_psi)
    
    logprec_w_psi <- logprec_w_psi_logrho_psi.update.details$logprec_w_psi
    logrho_psi <- logprec_w_psi_logrho_psi.update.details$logrho_psi
    cur.cormat_psi.inv <- logprec_w_psi_logrho_psi.update.details$cur.cormat_psi.inv
    cur.cormat_psi.logdet <- logprec_w_psi_logrho_psi.update.details$cur.cormat_psi.logdet
    cur.ss.w_psi <- logprec_w_psi_logrho_psi.update.details$cur.ss.w_psi
    att.logprec_w_psi_logrho_psi <- logprec_w_psi_logrho_psi.update.details$att.logprec_w_psi_logrho_psi
    acc.logprec_w_psi_logrho_psi <- logprec_w_psi_logrho_psi.update.details$acc.logprec_w_psi_logrho_psi
    
    prec_w_psi <-  exp(logprec_w_psi)
    s2_psi <- 1 / prec_w_psi
    s_psi <- sqrt(s2_psi)
    rho_psi <- exp(logrho_psi)
    
    # updated parameter: logprec_w_tau, logrho_tau
    
    logprec_w_tau_logrho_tau.update.details <- 
      logprec_w_tau_logrho_tau.update(nmesh, logprec_w_tau, logrho_tau, w_tau, 
                                      lambda_s_tau, lambda_rho_tau,
                                      cur.cormat_tau.inv, cur.cormat_tau.logdet,
                                      cur.ss.w_tau, c.mat, g1.mat, g2.mat,
                                      att.logprec_w_tau_logrho_tau, 
                                      acc.logprec_w_tau_logrho_tau, 
                                      mh.logprec_w_tau_logrho_tau)
    
    logprec_w_tau <- logprec_w_tau_logrho_tau.update.details$logprec_w_tau
    logrho_tau <- logprec_w_tau_logrho_tau.update.details$logrho_tau
    cur.cormat_tau.inv <- logprec_w_tau_logrho_tau.update.details$cur.cormat_tau.inv
    cur.cormat_tau.logdet <- logprec_w_tau_logrho_tau.update.details$cur.cormat_tau.logdet
    cur.ss.w_tau <- logprec_w_tau_logrho_tau.update.details$cur.ss.w_tau
    att.logprec_w_tau_logrho_tau <- logprec_w_tau_logrho_tau.update.details$att.logprec_w_tau_logrho_tau
    acc.logprec_w_tau_logrho_tau <- logprec_w_tau_logrho_tau.update.details$acc.logprec_w_tau_logrho_tau
    
    prec_w_tau <-  exp(logprec_w_tau)
    s2_tau <- 1 / prec_w_tau
    s_tau <- sqrt(s2_tau)
    rho_tau <- exp(logrho_tau)
    
    if(iter < (burn / 2)){
      this.update <- mhupdate(acc = acc.logprec_psi, att = att.logprec_psi, mh = mh.logprec_psi)
      acc.logprec_psi <- this.update$acc
      att.logprec_psi <- this.update$att
      mh.logprec_psi <- this.update$mh
      
      this.update <- mhupdate(acc = acc.logprec_tau, att = att.logprec_tau, mh = mh.logprec_tau)
      acc.logprec_tau <- this.update$acc
      att.logprec_tau <- this.update$att
      mh.logprec_tau <- this.update$mh
      
      this.update <- mhupdate(acc = acc.logprec_phi, att = att.logprec_phi, mh = mh.logprec_phi)
      acc.logprec_phi <- this.update$acc
      att.logprec_phi <- this.update$att
      mh.logprec_phi <- this.update$mh
      
      this.update <- mhupdate(acc = acc.logprec_w_psi_logrho_psi, 
                              att = att.logprec_w_psi_logrho_psi, 
                              mh = mh.logprec_w_psi_logrho_psi)
      acc.logprec_w_psi_logrho_psi <- this.update$acc
      att.logprec_w_psi_logrho_psi <- this.update$att
      mh.logprec_w_psi_logrho_psi <- this.update$mh
      
      this.update <- mhupdate(acc = acc.logprec_w_tau_logrho_tau, 
                              att = att.logprec_w_tau_logrho_tau, 
                              mh = mh.logprec_w_tau_logrho_tau)
      acc.logprec_w_tau_logrho_tau <- this.update$acc
      att.logprec_w_tau_logrho_tau <- this.update$att
      mh.logprec_w_tau_logrho_tau <- this.update$mh
    }
  }
    
    if(iter > burn){
      psi.sum <- psi.sum + psi
      psi2.sum <- psi2.sum + psi^2
      
      tau.sum <- tau.sum + tau
      tau2.sum <- tau2.sum + tau^2
      
      phi.sum <- phi.sum + phi
      phi2.sum <- phi2.sum + phi^2
      
      w_psi.sum <- w_psi.sum + w_psi
      w_psi2.sum <- w_psi2.sum + w_psi^2
      
      w_tau.sum <- w_tau.sum + w_tau
      w_tau2.sum <- w_tau2.sum + w_tau^2
    }
    
    # storage
    keepers.beta_psi[iter] <- beta_psi
    keepers.beta_tau[iter] <- beta_tau
    keepers.beta_phi[iter] <- beta_phi
    keepers.sigma_psi[iter] <- sigma_psi
    keepers.sigma_tau[iter] <- sigma_tau
    keepers.sigma_phi[iter] <- sigma_phi
    keepers.s_psi[iter] <- s_psi
    keepers.s_tau[iter] <- s_tau
    keepers.rho_psi[iter] <- rho_psi
    keepers.rho_tau[iter] <- rho_tau
    
    keepers.psi[iter, ] <- psi
    keepers.tau[iter, ] <- tau
    keepers.phi[iter, ] <- phi
    
    # if((iter %% 50 == 0)&(iter > 50)){
    #   par(mfrow = c(2, 5))
    #   plot(50:iter, keepers.beta_psi[50:iter], type = "l")
    #   plot(50:iter, keepers.beta_tau[50:iter], type = "l")
    #   plot(50:iter, keepers.beta_phi[50:iter], type = "l")
    #   plot(50:iter, keepers.sigma_psi[50:iter], type = "l")
    #   plot(50:iter, keepers.sigma_tau[50:iter], type = "l")
    #   plot(50:iter, keepers.sigma_phi[50:iter], type = "l")
    #   plot(50:iter, keepers.s_psi[50:iter], type = "l")
    #   plot(50:iter, keepers.s_tau[50:iter], type = "l")
    #   plot(50:iter, keepers.rho_psi[50:iter], type = "l")
    #   plot(50:iter, keepers.rho_tau[50:iter], type = "l")
    # }
    
    cat("\t iter", iter, "\n")
    
  } #end iters
  
  psi.posmean <- psi.sum / length(return.iters)
  psi.posvar <- psi2.sum / length(return.iters) - psi.posmean^2
  
  tau.posmean <- tau.sum / length(return.iters)
  tau.posvar <- tau2.sum / length(return.iters) - tau.posmean^2
  
  phi.posmean <- phi.sum / length(return.iters)
  phi.posvar <- phi2.sum / length(return.iters) - phi.posmean^2
  
  w_psi.posmean <- w_psi.sum / length(return.iters)
  w_psi.posvar <- w_psi2.sum / length(return.iters) - w_psi.posmean^2
  
  w_tau.posmean <- w_tau.sum / length(return.iters)
  w_tau.posvar <- w_tau2.sum / length(return.iters) - w_tau.posmean^2
  
  tock <- proc.time()[3]
  
  results <- list(beta_psi = keepers.beta_psi[return.iters],
                  beta_tau = keepers.beta_tau[return.iters],
                  beta_phi = keepers.beta_phi[return.iters],
                  sigma_psi = keepers.sigma_psi[return.iters], 
                  sigma_tau = keepers.sigma_tau[return.iters],
                  sigma_phi = keepers.sigma_phi[return.iters],
                  s_psi = keepers.s_psi[return.iters],
                  s_tau = keepers.s_tau[return.iters],
                  rho_psi = keepers.rho_psi[return.iters],
                  rho_tau = keepers.rho_tau[return.iters],
                  psi.selected = keepers.psi[return.iters, ],
                  tau.selected = keepers.tau[return.iters, ],
                  phi.selected = keepers.phi[return.iters, ],
                  psi.posmean = psi.posmean,
                  psi.posvar = psi.posvar,
                  tau.posmean = tau.posmean,
                  tau.posvar = tau.posvar,
                  phi.posmean = phi.posmean,
                  phi.posvar = phi.posvar,
                  w_psi.posmean = w_psi.posmean,
                  w_psi.posvar = w_psi.posvar,
                  w_tau.posmean = w_tau.posmean,
                  w_tau.posvar = w_tau.posvar,
                  minutes = (tock - tick) / 60
  )
  
  return(results)}
