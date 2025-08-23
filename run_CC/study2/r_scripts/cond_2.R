
rm(list=ls())

library(brms)   # for model fitting
library(MASS)   # for correlated predictors
source('fx.R')  # for simulation

# convenience function for diagnostics summary
diag_p = function(x) bayestestR::diagnostic_posterior(x, effects='all')

options(mc.cores=2) # 2 chains per model


# Specify constants (to be replaced) --------------------------------------
COND  = 2
CONST = c(15,5,1,0.09,0.3,0.09,0.3,0,0,0,0,0,0,0,0,0,1,1,1,0.4,0.4,0.4,0.1,0.9,100,2,3000,1000,1,0,0.64)
names(CONST) = c("N", "J", "G00", "TAU0", "G10", "TAU1", "G20", "TAU2", 
                 "G30", "TAU3", "RHO01", "RHO02", "RHO12", "X1_mu", 
                 "X2_mu", "X3_mu", "X1_sd", "X2_sd", "X3_sd",
                 "R_12","R_13","R_23", "TAUX","SIGMAX",
                 "S","nchain","niter","burn","thin","marg","SIG2")
CONST = as.data.frame(t(CONST))

  
# Compile models ----------------------------------------------------
# Fit initial model to dummy data, then just call update below to avoid recompiling time.
# This is OS-specific, so has to be run once on the cluster before each simulation cell.
# Fit five models:
# A. True model
# B. The (incorrect) absence of a random effect (TAU1 = 0)
# C. The (incorrect) presence of a random effect (TAU2 > 0) 
# D. The (incorrect) absence of a fixed effect (G20 = 0)
# E. The (incorrect) presence of a fixed effect (G30 != 0)

DM    = expand.grid(1:10,1:10) # design matrix (dummy N & J)
DM$X1 = rnorm(nrow(DM)) # mean and sd don't matter here, just setting up model structure
DM$X2 = rnorm(nrow(DM))
DM$X3 = rnorm(nrow(DM))
DM$y  = rnorm(nrow(DM)) # same for the outcome
colnames(DM) = c('obs','g','X1','X2','X3', 'y')

options(mc.cores=CONST$nchain)

modA = brm(y ~ 1 + X1 + X2 + (X1|g),      data=DM, chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin, silent = 2)
modB = brm(y ~ 1 + X1 + X2 + (1|g),       data=DM, chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin, silent = 2)
modC = brm(y ~ 1 + X1 + X2 + (X1+X2|g),   data=DM, chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin, silent = 2)
modD = brm(y ~ 1 + X1 + (X1|g),           data=DM, chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin, silent = 2)
modE = brm(y ~ 1 + X1 + X2 + X3 + (X1|g), data=DM, chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin, silent = 2)



# Fit each model to each dataset and compute fit indices --------------------------------------------------

out   = list() # output to save

for(s in 1:CONST$S) {
  
  cat('\n=========================generating data for sample ', s,'=========================\n')
  
  # Simulate from true model each sample
  g  = rep(1:CONST$N, each=CONST$J)
  
  ## Produce correlated predictors with particular ICC 
  MU  = c(CONST$X1_mu, CONST$X2_mu, CONST$X3_mu)
  L1SIG = matrix(c(CONST$SIGMAX, CONST$R_12*sqrt(CONST$SIGMAX)*sqrt(CONST$SIGMAX), CONST$R_13*sqrt(CONST$SIGMAX)*sqrt(CONST$SIGMAX),
                   CONST$R_12*sqrt(CONST$SIGMAX)*sqrt(CONST$SIGMAX), CONST$SIGMAX, CONST$R_23*sqrt(CONST$SIGMAX)*sqrt(CONST$SIGMAX),
                   CONST$R_13*sqrt(CONST$SIGMAX)*sqrt(CONST$SIGMAX), CONST$R_23*sqrt(CONST$SIGMAX)*sqrt(CONST$SIGMAX), CONST$SIGMAX), 
                 3, 3, byrow=T)
  L2SIG =  matrix(c(CONST$TAUX, CONST$R_12*sqrt(CONST$TAUX)*sqrt(CONST$TAUX), CONST$R_13*sqrt(CONST$TAUX)*sqrt(CONST$TAUX),
                    CONST$R_12*sqrt(CONST$TAUX)*sqrt(CONST$TAUX), CONST$TAUX, CONST$R_23*sqrt(CONST$TAUX)*sqrt(CONST$TAUX),
                    CONST$R_13*sqrt(CONST$TAUX)*sqrt(CONST$TAUX), CONST$R_23*sqrt(CONST$TAUX)*sqrt(CONST$TAUX), CONST$TAUX), 
                  3, 3, byrow=T)
  
  
  Xl1 = rmvnorm(CONST$N*CONST$J, mean=MU, sigma = L1SIG) # l1
  Xl2 = rmvnorm(CONST$N, mean=MU, sigma = L2SIG)         # l2
  Xl2_ = lapply(1:nrow(Xl2), function(i){matrix(rep(Xl2[i,], CONST$J), ncol=3, byrow=TRUE)})
  Xl2 = do.call("rbind",Xl2_)
  
  dat = data.frame(X1 = Xl2[,1]+Xl1[,1], 
                   X2 = Xl2[,2]+Xl1[,2], 
                   X3 = Xl2[,3]+Xl1[,3], 
                   g  = g)
  
  # Store pars in list
  pars = list(GAMMA=c(CONST$G00, CONST$G10, CONST$G20), 
              TAU=matrix(c(CONST$TAU0, 
                           CONST$RHO01*sqrt(CONST$TAU0)*sqrt(CONST$TAU1), 
                           CONST$RHO01*sqrt(CONST$TAU0)*sqrt(CONST$TAU1), 
                           CONST$TAU1), 2, 2), 
              SIGMA=CONST$SIG2)
  names(pars$GAMMA) = c('(Intercept)','X1','X2')
  
  # simulate from true model S times
  dat$y = simulate_lmer(y ~ X1 + X2 + (X1|g), DM=dat, pars=pars, nsim = 1, seed = runif(1,0,10000))[,1]
  
  cat('\n=========================fitting sample ', s,'=========================\n')
  
  # fit
  cat('fitting models...\n')
  suppressMessages(
    {
      modA = update(modA, newdata=dat, recompile = F, refresh=0)
      modB = update(modB, newdata=dat, recompile = F, refresh=0)
      modC = update(modC, newdata=dat, recompile = F, refresh=0)
      modD = update(modD, newdata=dat, recompile = F, refresh=0)
      modE = update(modE, newdata=dat, recompile = F, refresh=0)
    }
  )
  
  # Compute fit indices for each model 
  cat('computing fit indices...\n')
  
  suppressWarnings(
    {
      dicA  = DIC_brms(modA)
      waicA = WAIC(modA)$estimates
      looA  = LOO(modA)$estimates
      
      dicB  = DIC_brms(modB)
      waicB = WAIC(modB)$estimates
      looB  = LOO(modB)$estimates
      
      dicC  = DIC_brms(modC)
      waicC = WAIC(modC)$estimates
      looC  = LOO(modC)$estimates
      
      dicD  = DIC_brms(modD)
      waicD = WAIC(modD)$estimates
      looD  = LOO(modD)$estimates
      
      dicE  = DIC_brms(modE)
      waicE = WAIC(modE)$estimates
      looE  = LOO(modE)$estimates
    }
  )
  
  # Save output
  # (saving the whole models would take too much memory)
  # (same with the pll for now)
  cat('storing output...\n')
  
  out[[s]] = list(cond=COND, iter=s, marg=CONST$marg,
             data  = dat,
             # pll = list(A=log_lik(modA),B=log_lik(modB),C=log_lik(modC),D=log_lik(modD),E=log_lik(modE)),
             fixef = list(A=fixef(modA),B=fixef(modB),C=fixef(modC),D=fixef(modD),E=fixef(modE)),
             ranef = list(A=VarCorr(modA),B=VarCorr(modB),C=VarCorr(modC),D=VarCorr(modD),E=VarCorr(modE)),
             diag  = list(A=diag_p(modA),B=diag_p(modB),C=diag_p(modC),D=diag_p(modD),E=diag_p(modE)),
             dic   = list(A=dicA,B=dicB,C=dicC,D=dicD,E=dicE),
             waic  = list(A=waicA,B=waicB,C=waicC,D=waicD,E=waicE), 
             loo   = list(A=looA,B=looB,C=looC,D=looD,E=looE)
             )
  
} # end of sample loop 

# save output to disk again at end
fname = paste0('/home/sdevine/scratch/cond_',COND,'.rds')
saveRDS(out, fname)

# Sys.sleep(500) # give ample time to save to avoid data loss


