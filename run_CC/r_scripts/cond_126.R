
rm(list=ls())

library(brms)   # for model fitting
source('fx.R')  # for simulation

# convenience function for diagnostics summary
diag_p = function(x) bayestestR::diagnostic_posterior(x, effects='all')

options(mc.cores=2) # 2 chains per model


# Specify constants (to be replaced) --------------------------------------
COND  = 126
CONST = c(70,100,1,0.16,0.2,0.16,0.4,0,0,0,0.1,0,0,0,0,0,1,1,1,100,2,3000,1000,1,0.48)
names(CONST) = c("N", "J", "G00", "TAU0", "G10", "TAU1", "G20", "TAU2", 
                 "G30", "TAU3", "RHO01", "RHO02", "RHO12", "X1_mu", 
                 "X2_mu", "X3_mu", "X1_sd", "X2_sd", "X3_sd", "S", 
                 "nchain","niter","burn","thin", "SIG2")
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


# Generate dataset --------------------------------------------------------

# Simulate from true model
g  = rep(1:CONST$N, each=CONST$J)
X1 = rnorm(n=CONST$N*CONST$J, mean=CONST$X1_mu, sd=CONST$X1_sd)
X2 = rnorm(n=CONST$N*CONST$J, mean=CONST$X2_mu, sd=CONST$X2_sd)
X3 = rnorm(n=CONST$N*CONST$J, mean=CONST$X3_mu, sd=CONST$X3_sd)
dm = data.frame(X1=X1,X2=X2,X3=X3,g=g)

# Store pars in list
pars = list(GAMMA=c(CONST$G00, CONST$G10, CONST$G20), 
            TAU=matrix(c(CONST$TAU0, 
                         CONST$RHO01*sqrt(CONST$TAU0)*sqrt(CONST$TAU1), 
                         CONST$RHO01*sqrt(CONST$TAU0)*sqrt(CONST$TAU1), 
                         CONST$TAU1), 2, 2), 
            SIGMA=CONST$SIG)
names(pars$GAMMA) = c('(Intercept)','X1','X2')

# simulate from true model S times
y = simulate_lmer(y ~ X1 + X2 + (X1|g), DM=dm, pars=pars, nsim = CONST$S)

# Fit each model to each dataset and compute fit indices --------------------------------------------------

out   = list() # output to save

for(s in 1:CONST$S) {
  
  cat('\n=========================fitting sample ', s,'=========================\n')
  
  dat = dm
  dat$y = y[,s]
  
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
  
  out[[s]] = list(cond=COND, iter=s, 
             data  = dm,
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


