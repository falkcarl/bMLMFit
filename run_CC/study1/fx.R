
######################################################################
##
## Copyright 2024 Carl F. Falk and Sean Devine
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
## <http://www.gnu.org/licenses/>


library(mnormt) # for simulations
library(MCMCpack)
library(matrixcalc)
library(brms) 
library(lme4)
library(dplyr)
library(Matrix)
library(numDeriv)
library(mvtnorm)
library(loo)

#' Simulate one lmer response from a formula and design matrix
#' @param formula: R formula object, in lme4 style
#' @param DM: Design matrix, with all relevant variables present
#' @param seed: Random seed.
#'  
#' @returns vector with simulated outcome
sim1 = function(formula, DM, pars, seed=2022) {

  set.seed(seed)  
  
  # 1. Extract pars
  gam = pars$GAMMA
  v   = pars$TAU
  sig = pars$SIGMA
  
  # 2. Extract vars from formula
  preds = attributes(terms(formula))$term.labels
  f     = preds[-length(preds)]
  r     = gsub(' ', '', strsplit(preds[length(preds)],'\\|')[[1]])
  g     = r[length(r)]
  r     = r[-length(r)]
  
  dmat = cbind(DM[,g], 1, DM[,f])
  n    = length(unique(dmat[,1]))
  
  # 3. Sample random effects (level-2)
  Uj = rmnorm(n, mean = rep(0, length(r)+1), varcov = v) # intercept, slopes
  
  # 4. Sample residuals (level-1)
  Rij = rnorm(nrow(dmat), 0, sqrt(sig))
  
  # 5. Compute betas as GAMMA + Uj
  # First, add zero column for fixed effects only
  gam_m = matrix(gam, n,length(gam), byrow = TRUE)
  if(ncol(Uj)!=ncol(gam_m)) {
    diff = ncol(gam_m) - ncol(Uj)
    if(diff < 0) stop('VarCov mispecified')
    for(i in 1:diff) Uj = cbind(Uj, 0)
  }
  b = gam_m + Uj
  
  # 6. Compute response
  y = rowSums(dmat[,-1]*b[dmat[,1],]) + Rij
  
  y
  
}

simulate_lmer = function(formula, DM, pars, nsim, seed=2022) {
  #' Wrapper for sim1(), which reproduces the process `nsim` times
  
  set.seed(seed)
  sapply(1:nsim, function(i) sim1(formula, DM, pars, runif(1,0,1e6)))
  
}

#' Convenience function to compute the DIC for a brms object
#' Computed at the mean posterior parameter estimates, as implemented in Gelman et al., 2014, Stats. & Comp.
#' @param mod: a fitted brms object
#' @param version: how to compute pDIC (see Gelman et al., 2014, Stats. & Comp.)
#' @return [vector] expected log predictive density, effective number of parameters, dic value
DIC_brms = function(mod, version=1) {
  

  # p(y | MAP)
  y    = mod$data[,1]
  yhat = predict(mod)[,'Estimate']
  sig  = summary(mod)$spec_pars[,'Estimate']
  pMAP = sum(dnorm(y, yhat, sig, log=T))
  
  # pDIC
  ll = log_lik(mod) # point-wise likelihood
  
  if(version==1) {
    mll   = mean(rowSums(ll))    # mean deviance
    pDIC  = 2*(pMAP-mll)
  } else {
    pDIC = 2*var(rowSums(ll))
  }
  
  # elpd
  elpd = pMAP - pDIC
  
  # compute dic
  return(c(elpd = elpd, pDIC=pDIC, DIC=-2*pMAP + 2*pDIC))
  
}


#' Computes error variance, such that total variance sums to 1
#' Note: this is specific for these simulations, so only takes arguments
#' necessary. I.e., it is not generalizable in its current form.
#'
#' @param g00:    fixed intercept
#' @param g10:    fixed slope 1
#' @param g20:    fixed slope 2
#' @param varx1:  variance of X1
#' @param varx2:  variance of X2
#' @param cov12:  covariance in X1 and X2
#' @param tau0:   random intercept
#' @param tau1:   random slope for X1
#' @param cov01:  covariance between tau0 and tau1
#' @param return_prop_var: Return proportions of variances explained
#'
#'@return error variance
get_vareY = function (gam0,gam1,gam2,varx1,varx2,cov12,tau0,tau1,cov01, 
                      return_prop_var=F) {
  
  # make random effects matrix
  ranmat = makeran(
    ty    = tau0,
    tgam1 = tau1,
    covintgam1 = cov01*sqrt(tau0)*sqrt(tau1)
  )
  
  # tau00 (Eq 14)
  varintY = ranmat["ty","ty"]
  
  # T - covariance among random effects, not including intercept
  #vargam12<-ranmat[c("tgam1","tgam2"),c("tgam1","tgam2")]
  
  ###################################
  # computations: Y
  
  # within-group variance of X
  # Phi^w
  covX = matrix(c(varx1,cov12,cov12,varx2),2,2,byrow=T)
  covW = cbind(c(0,0,0),rbind(c(0,0),covX))
  
  # variance due to fixed level 1 predictors (Eq. 11)
  varfixed = c(gam1, gam2) %*%covX %*% c(gam1, gam2)
  
  # variance due to level 1 predictors via random slope variation (Eq. 13)
  #varrand<- sum(diag(vargam12%*%covX)) # actually reduces to this due to 0's in covW (i.e., Phi^w)
  varrand = sum(diag(ranmat%*%covW))
  
  # Note: vareY is Eq. 15, variance of level 1 residuals
  
  # Total variance of outcome (Eq. 10)
  
  # return error variance so that total variance sums to 1
  vareY = 1-(varintY+varfixed+varrand)
  
  if(return_prop_var) {
    vary<-varintY+varfixed+varrand+vareY
    
    return(list(Y.fixed = as.numeric(varfixed/vary),
              Y.rand    = as.numeric(varrand/vary),
              Y.int     = as.numeric(varintY/vary),
              Y.e       = as.numeric(vareY/vary),
              varY=vary)
    )
  }
  
  vareY[[1]]
  
}


# Creates variance-covariance matrix of random effects for 2 predictor regression model (matrix that contains tau's)
# ty - random intercept variance
# tgam1 - random variance of first slope
# tgam2 - random variance of second slope (currently 0 by default)
# covintgam1 - covariance between random intercept and random first slope (Sean has this labeled rho01)
# covintgam2 - covariance between random intercept and random first slope (currently zero)
# covgam1gam2 - covariance between random slopes (currently 0)
makeran<-function(ty=.25,tgam1=.25,tgam2=0,covgam1gam2=0,covintgam1=.1,covintgam2=0){

  ran<-matrix(c(
    ty, covintgam1, covintgam2, #dy
    covintgam1, tgam1, covgam1gam2, #gam1
    covintgam2, covgam1gam2, tgam2 # gam2
  ),nrow=3,ncol=3,byrow=TRUE)
  colnames(ran)<-rownames(ran)<-c("ty","tgam1","tgam2")
  stopifnot(isSymmetric(ran))
  stopifnot(is.positive.semi.definite(ran))
  return(ran)
}


# Set up model matrices for MLM - for use with custom marginal likelihood code
# in R or stan
#
# Note: written to be kind of general, but not completely general for all lmer-style syntax
# Does not allow for || in random effects (i.e., uncorrelated random effects)
# Assumes model has an intercept and random intercept.
# Must have fixed effect for us to have it in there as a random effect, and so on
# Only really tested with models for a simulation study.
mod.mat.mlm <- function(formula, data, sparse=FALSE){
  
  n <- nrow(data) # sample size
  In <- diag(n) # identity matrix of size n x n
  
  # helper function from lme4 to extract stuff about random effects
  re <- lme4::findbars(formula)
  clustvar <- all.vars(re[[1]][[3]])
  zs <- all.vars(re[[1]][[2]])
  
  # ensure data is sorted so that random effects matrix does not break
  if(is.unsorted(data[,clustvar])){
    stop("TODO: fix what happens if data is not sorted by cluster variable")
  }
  clustsizes <- table(data[,clustvar]) # cluster sizes; useful for stan
  
  # y
  y <- deparse(formula[[2]])
  y <- data[,y]
  
  # Create design matrices
  
  # X - for fixed effects
  xs <- all.vars(formula[[c(3,2)]])
  Xmat <- cbind(1, data[,xs])
  
  # Z - for random effects
  # Design matrix for random effects is a bit tricky to put together
  # We want a block-diagonal matrix; each block has a little design matrix for each participant
  Zmat <- bdiag(lapply(unique(data[,clustvar]), function(x){
    datsub <- data[data[,clustvar]==x,]
    if(length(zs)>0){
      out <- as.matrix(cbind(1,datsub[,zs]))
    } else {
      out <- matrix(1,ncol = 1, nrow = nrow(datsub))
    }
    out
  }))
  # if we don't actually want a sparse matrix
  if(!sparse){
    Zmat <- as.matrix(Zmat)
  }
  
  nfixef <- 1 + length(xs) # number of variables + intercept
  p <- 1+length(zs)
  nranef <- (p)*(p+1)/2 # unique elements in random effects matrix,  p*(p+1)/2 for p-dimensional matrix
  
  return(list(n = n,
              nclust = length(unique(data[,clustvar])),
              nfixef = nfixef,
              nranef = nranef,
              p = p,
              y = y,
              In = In,
              Xmat = as.matrix(Xmat),
              Zmat = Zmat,
              #xs = xs, # could be useful later if we want pretty output; i.e., match fixed and random effects to text names for their coefficients
              #zs = zs,
              clustsizes = clustsizes))
}
# test
#mats <- mod.mat.mlm(y ~ X1 + X2 + (1|g), data=dat)


# log-likelihood for MLM
# See readings from Bates et al and Wang & Merkle (2018) for the log-likelihood equation
#
# par - parameters
# y - outcome
# X and Z - design matrices
# In - identity matrix
# clustid - level 2 indices; not exactly used, but I needed to know number of level 2 units (i.e., number of clusters)
ll <- function(par, # vector of parameters
               nfixef, nranef, # number of each type of parameter
               y, X, Z, In, # model matrices
               nclust, # number of clusters
               neg=TRUE){
  n <- length(y)
  sig2 <- par[1]    # error variance    
  betas <- par[2:(nfixef+1)] # fixed effect coefficients
  gs <- par[(nfixef+2):(nfixef+nranef+1)] # random effect variances
  
  XB <- X %*% betas
  R <- sig2*In
  G <- bdiag(lapply(1:nclust, function(x){xpnd(gs)}))
  
  V <- Z%*%G%*%t(Z) + R
  
  # solve(V) # also trying cholesky decomp instead in case it's more numerically stable
  out <- -(.5*n)*log(2*pi) - (.5)*determinant(V)$modulus - .5*t(y-XB)%*%chol2inv(chol(V))%*%(y-XB)
  if(neg){out <- -out}
  
  return(as.numeric(out))
}


# Re-do of estimation using marginal likelihood
# Wrapper function
fit.mlm.marg <- function(data, formula, estimator = "nlminb"){
  
  mats <- mod.mat.mlm(formula, data=dat)
  
  # starting values are still a problem...  
  par.start <- c(var(mats$y), mean(mats$y),
                 rep(0, mats$nfixef-1),
                 vech(diag(rep(.05,mats$p))))
  
  #par.start <- c(.7, 1, .1, .4, .03, .01, .03)
  
  # solve
  if(estimator == "nlm"){
    
    nlm.res <- nlm(ll, par.start,
                   nfixef = mats$nfixef, nranef = mats$nranef,
                   y=mats$y, X=mats$Xmat, Z=mats$Zmat, In=mats$In,
                   nclust = mats$nclust,
                   hessian=TRUE)
    
    est <- nlm.res$estimate
    ses <- sqrt(diag(solve(nlm.res$hessian)))
    
    # arrange results
    out <- list(est = est,
                ses = ses,
                code = nlm.res$code,
                it = nlm.res$iterations,
                nll = nlm.res$minimum)
    
    
  } else if (estimator == "nlminb"){
    
    nlb.res <- nlminb(par.start, ll, nfixef = mats$nfixef, nranef = mats$nranef,
                      y=mats$y, X=mats$Xmat, Z=mats$Zmat, In=mats$In,
                      nclust = mats$nclust)
    
    H<-hessian(ll, nlb.res$par, nfixef = mats$nfixef, nranef = mats$nranef,
               y=mats$y, X=mats$Xmat, Z=mats$Zmat, In=mats$In,
               nclust = mats$nclust)
    
    est <- nlb.res$par
    ses <- sqrt(diag(solve(H)))
    
    # arrange results
    out <- list(est = est,
                ses = ses,
                code = nlb.res$convergence,
                it = nlb.res$iterations,
                nll = nlb.res$objective)
  }
  
  return(out)
  
}


#' CFF modified from Sean's code here: https://github.com/seandamiandevine/bMLMFit
#' Convenience function to compute the DIC for a brms object
#' Computed at the mean posterior parameter estimates, as implemented in Gelman et al., 2014, Stats. & Comp.
#' @param mod: a fitted stan object, assuming mlmmarg.stan was used
#' @param stan.dat: model matrices, e.g., using mod.mat.mlm, with addition of cluster sizes
#' @param version: how to compute pDIC (see Gelman et al., 2014, Stats. & Comp.)
#' @return [vector] expected log predictive density, effective number of parameters, dic value
DIC_marg = function(mod, stan.dat, version=1) {

  # p(y | MAP)
  y    = stan.dat$y
  betas <- (summary(mod, pars="betas")$summary)[,"mean"] # grab beta point estimates
  yhat = stan.dat$Xmat %*% betas # prediction with model matrices
  sig2  =  (summary(mod, pars="sig")$summary)[,"mean"]^2 # grab sigma point estimate from model
  G <- (summary(mod, pars="Tau")$summary)[,"mean"] # grab random effect covariance matrix here
  #G <- xpnd(G) #Is this needed, or just coerce to a matrix? Try w/ 2 random effects to test
  G <- matrix(G, nrow=stan.dat$p, ncol=stan.dat$p)
  
  pMAP <- 0
  idx<-1
  idx2<-1
  for(n in stan.dat$clustsizes){
    Z <- stan.dat$Zmat[idx:(idx+n-1),idx2:(idx2+stan.dat$p-1)]
    R <- diag(rep(sig2, n))
    pMAP <- pMAP + dmvnorm(y[idx:(idx+n-1)], yhat[idx:(idx+n-1)], Z%*%G%*%t(Z) + R, log=TRUE)
    idx <- idx + n
    idx2 <- idx2 + stan.dat$p
  }
  
  # pDIC
  ll = extract_log_lik(mod) # point-wise likelihood
  
  if(version==1) {
    mll   = mean(rowSums(ll))    # mean deviance
    pDIC  = 2*(pMAP-mll)
  } else {
    pDIC = 2*var(rowSums(ll))
  }
  
  # elpd
  elpd = pMAP - pDIC
  
  # compute dic
  return(c(elpd = elpd, pDIC=pDIC, DIC=-2*pMAP + 2*pDIC))
  
}

#' Function to summarize parameters of interest from stanfit
#' @param mod a stanfit object
#' @param parname name of parameter to extract
summarize_post = function(mod,parname) {
  
  post_sum = summary(mod)$summary
  post_sum = post_sum[grepl(parname,rownames(post_sum)),c('mean','se_mean','2.5%','97.5%')]
  
  if(is.null(ncol(post_sum))) post_sum = matrix(post_sum, 1, length(post_sum))
  
  colnames(post_sum) = c("Estimates","Est.Error","Q2.5","Q97.5")

  return(post_sum)
}

#' Convenience function that mimics brms' fixef() for a MLM fit with marginal estimation
#' @param mod a stanfit object
fixef_marg = function(mod) {
  post_sum = summarize_post(mod, 'betas')
  rownames(post_sum) = c("Intercept",paste0("X",1:(nrow(post_sum)-1)))
  post_sum
}

#' Convenience function that mimics brms' VarCorr() for a MLM fit with marginal estimation
#' @param mod a stanfit object
ranef_marg = function(mod) {
  
  out = list()
  
  out[['g']] = list(sd=summarize_post(mod,'sig'),
                  tau=summarize_post(mod,'Tau'), 
                  cor=summarize_post(mod,'gsmat'))
  
  out
  
}

#' Convenience function to extract diagnostic measures for a MLM fit with marginal estimation
#' @param mod a stanfit object
diag_marg = function(mod) {
  
  post_sum = summary(mod)$summary
  post_sum = post_sum[!grepl('yhat|log_lik|idx|idxzmat|lp__',rownames(post_sum)),c('n_eff','Rhat')]
  colnames(post_sum) = c('ESS','Rhat')
  
  post_sum  
  
}


