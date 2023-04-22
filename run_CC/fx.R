
library(mnormt) # for simulations
library(MCMCpack)
library(matrixcalc)
library(brms) 

sim1 = function(formula, DM, pars, seed=2022) {
  #' Simulate one lmer response from a formula and design matrix
  #' @param formula: R formula object, in lme4 style
  #' @param DM: Design matrix, with all relevant variables present
  #' @param seed: Random seed.
  #'  
  #' @returns vector with simulated outcome
  
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

DIC_brms = function(mod, version=1) {
  
  #' Convenience function to compute the DIC for a brms object
  #' Computed at the mean posterior parameter estimates, as implemented in Gelman et al., 2014, Stats. & Comp.
  #' @param mod: a fitted brms object
  #' @param version: how to compute pDIC (see Gelman et al., 2014, Stats. & Comp.)
  #' @return [vector] expected log predictive density, effective number of parameters, dic value
  
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


get_vareY = function (gam0,gam1,gam2,varx1,varx2,cov12,tau0,tau1,cov01, 
                      return_prop_var=F) {
  
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



makeran<-function(ty=.25,tgam1=.25,tgam2=0,covgam1gam2=0,covintgam1=.1,covintgam2=0){
  # Creates variance-covariance matrix of random effects for 2 predictor regression model (matrix that contains tau's)
  # ty - random intercept variance
  # tgam1 - random variance of first slope
  # tgam2 - random variance of second slope (currently 0 by default)
  # covintgam1 - covariance between random intercept and random first slope (Sean has this labeled rho01)
  # covintgam2 - covariance between random intercept and random first slope (currently zero)
  # covgam1gam2 - covariance between random slopes (currently 0)
  
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

