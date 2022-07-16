rm(list = ls())

mydir <- "C:/Users/Arnab Hazra/Desktop/max_and_smooth"
setwd(mydir)

x <- y <- seq(0.25, 5, 0.25)

loc <- as.matrix(expand.grid(x=x, y=y))
colnames(loc) <- c("lon", "lat")

loc.scaled <- apply(loc, 2, scale)

mu <- 10 + 0.5 * loc.scaled[ , 1] - 0.5 * loc.scaled[ , 2]
sigma <- exp(1 + 0.5 * loc.scaled[ , 1]^2 - 0.5 * loc.scaled[ , 2]^2)

ns <- nrow(loc)
nt <- 1000 # number of temporal replications

Y <- mu + sigma * matrix(rnorm(ns * nt), ns, nt) # simulated data matrix

dim(Y)

#------------------
# Stage 1

source("stage1_functions.R")

thresholds <- apply(Y, 1, function(x){quantile(x[x > 0], probs = 0.90)})

ppfit.mles <- lapply(1:ns, ppfit.mle, thresholds = quantiles90, Y = Y)

mles <- t(sapply(seq_along(ppfit.mles), function(station.no){ppfit.mles[[station.no]]$mle}))

# next, we transform to the original scale

mles.original <- mles

mles.original[ , 1] <- exp(mles[ , 1])
mles.original[ , 2] <- exp(mles[ , 1] + mles[ , 2])
mles.original[ , 3] <- g(mles[ , 3])
npy <- 365.25 # number of days per year

covmats <- lapply(1:ns, rPPP) # this step is time consuming

mles.covmats <- lapply(1:ns, function(station.no){
  list(mle = ppfit.mles[[station.no]]$mle, covmat = covmats[[station.no]])})

# store the outputs from Stage 1

save(mles.covmats, file = "mles_covmats_final.Rdata")

#------------------
# Stage 2

source("stage2_functions.R")

library(INLA)

mesh.domain <- inla.mesh.2d(loc = loc, max.n.strict = 250, offset = 0.08, 
                            max.edge = c(0.1, 1), cutoff = 0.8)

A <- inla.spde.make.A(mesh = mesh.domain, loc = loc)
fem.mesh <- inla.mesh.fem(mesh.domain, order = 2)
c.mat <- fem.mesh$c0
g1.mat <- fem.mesh$g1
g2.mat <- fem.mesh$g2
inla.mats <- list(c.mat = c.mat, g1.mat = g1.mat, g2.mat = g2.mat, A = A)

plot(mesh.domain)

fit.mcmc.maxNsmooth <- mcmc.maxNsmooth(mles.covmats, loc, inla.mats, alpha = 2,
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
                                       iters = 5000, burn = 1000, thin = 2)

# store the outputs from Stage 2

save(fit.mcmc.maxNsmooth, file = "fit.mcmc.maxNsmooth.Rdata")

#------------------
# return level calculation

mu.chain <- exp(fit.mcmc.maxNsmooth$psi.selected)
sigma.chain <- exp(fit.mcmc.maxNsmooth$psi.selected + fit.mcmc.maxNsmooth$tau.selected)
xi.chain <- g(fit.mcmc.maxNsmooth$phi.selected)

# M-year return level

M <- 20

rl.chain <- mu.chain + sigma.chain / xi.chain * ((-log(1 - 1 / M))^(-xi.chain) - 1)

rl.posmean <- apply(rl.chain, 2, mean)
rl.possd <- apply(rl.chain, 2, sd)

library(plot3D)

par(mfrow = c(1, 2))
scatter2D(x = loc[ , 1], y = loc[ , 2], colvar = rl.posmean, pch = 15, main = "Posterior mean")
scatter2D(x = loc[ , 1], y = loc[ , 2], colvar = rl.possd, pch = 15, main = "Posterior SD")
