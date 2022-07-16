details of the files:

#------------------
stage1_functions.R
#------------------

g, h: transformation functions

pp.fit.fast: estimates the point process parameters (transformed, 
	     log-mu, log-sigma/mu, phi) for each grid cell, also
	     returns the covariance matrices (negative Hessian) 
	     of the estimators, we keep only the estimates later

ppfit.mle: calls pp.fit.fast to perform for each data location

rPPP: using 1000 parametric bootstrap samples, calculates the covariance matrices of the estimators

mles_covmats_final.Rdata: mles and their covariance matrices across locations obtained from the 
			  outputs of functions ppfit.mle and rPPP

#------------------
stage2_functions.R
#------------------

mhupdate: adaptive M-H algorithm

cormat.update: calculates dense Matern covariance matrix

cormat.inv.update.inla: calculates sparse Matern covariance matrix

eta.update: updates eta

nu.update: updates nu

logprec_psi.update: updates sigmaSq_psi (f(x) = log(1/x) transformed)

logprec_tau.update: updates sigmaSq_tau (f(x) = log(1/x) transformed)

logprec_phi.update: updates sigmaSq_phi (f(x) = log(1/x) transformed)

logprec_w_psi_logrho_psi.update: updates sSq_psi (f(x) = log(1/x) transformed) 
				 and rho_psi (f(x) = log(x) transformed) jointly

logprec_w_tau_logrho_tau.update: updates sSq_tau (f(x) = log(1/x) transformed) 
				 and rho_tau (f(x) = log(x) transformed) jointly

mcmc.maxNsmooth: performs MCMC

#------------------
example.R
#------------------

Here we create an artifical dataset on a 40X40 regular grid from independent normal 
distributions. call it Y, and the set of data locations are in the 1600X2 matrix loc

Once the users have a spatiotemporal data matrix (denote by Y) and coordinates (denoted by loc),
they can simply follow stage 1 and stage 2

Based on the MCMC chains for log-mu, log-sigma/mu, phi, return level maps can be easily 
obtained after transformation




